/*
 * Copyright (C) 2010-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "VirtualSitesRelative.hpp"

#ifdef VIRTUAL_SITES_RELATIVE

#include "Particle.hpp"
#include "cells.hpp"
#include "forces.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "rotation.hpp"

#include <utils/Vector.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/math/tensor_product.hpp>
#include <utils/quaternion.hpp>

#include <stdexcept>

namespace {
/**
 * @brief Vector pointing from the real particle to the virtual site.
 *
 * @return Relative distance.
 */
auto connection_vector(
    Particle const &p_ref,
    ParticleProperties::VirtualSitesRelativeParameters const &vs_rel) {
  // Calculate the quaternion defining the orientation of the vector connecting
  // the virtual site and the real particle
  // This is obtained, by multiplying the quaternion representing the director
  // of the real particle with the quaternion of the virtual particle, which
  // specifies the relative orientation.
  auto const director = Utils::convert_quaternion_to_director(
                            p_ref.quat() * vs_rel.rel_orientation)
                            .normalize();

  return vs_rel.distance * director;
}

/**
 * @brief Velocity of the virtual site
 * @param p_ref Reference particle for the virtual site.
 * @param vs_rel Parameters for the virtual site.
 * @return Velocity of the virtual site.
 */
Utils::Vector3d
velocity(Particle const &p_ref,
         ParticleProperties::VirtualSitesRelativeParameters const &vs_rel) {
  auto const d = connection_vector(p_ref, vs_rel);

  // Get omega of real particle in space-fixed frame
  auto const omega_space_frame =
      convert_vector_body_to_space(p_ref, p_ref.omega());
  // Obtain velocity from v=v_real particle + omega_real_particle \times
  // director
  return vector_product(omega_space_frame, d) + p_ref.v();
}

/**
 * @brief Get reference particle.
 *
 * @param vs_rel Parameters to get the reference particle for.
 * @return Pointer to reference particle.
 */
Particle &get_reference_particle(
    ParticleProperties::VirtualSitesRelativeParameters const &vs_rel) {
  auto p_ref_ptr = cell_structure.get_local_particle(vs_rel.to_particle_id);
  if (!p_ref_ptr) {
    throw std::runtime_error("No real particle associated with virtual site.");
  }
  return *p_ref_ptr;
}

/**
 * @brief Constraint force to hold the particles at its prescribed position.
 *
 * @param f Force on the virtual site.
 * @param p_ref Reference particle.
 * @param vs_rel Parameters.
 * @return Constraint force.
 */
auto constraint_stress(
    const Utils::Vector3d &f, const Particle &p_ref,
    const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  /* The constraint force is minus the force on the particle, make it force
   * free. The counter force is translated by the connection vector to the
   * real particle, hence the virial stress is */
  return tensor_product(-f, connection_vector(p_ref, vs_rel));
}
} // namespace

void VirtualSitesRelative::update() const {
  cell_structure.ghosts_update(Cells::DATA_PART_POSITION |
                               Cells::DATA_PART_MOMENTUM);

  auto const particles = cell_structure.local_particles();
  for (auto &p : particles) {
    if (!p.is_virtual())
      continue;

    auto const &p_ref = get_reference_particle(p.vs_relative());

    auto new_pos = p_ref.pos() + connection_vector(p_ref, p.vs_relative());
    /* The shift has to respect periodic boundaries: if the reference
     * particles is not in the same image box, we potentially avoid shifting
     * to the other side of the box. */
    auto shift = box_geo.get_mi_vector(new_pos, p.pos());
    p.pos() += shift;
    Utils::Vector3i image_shift{};
    fold_position(shift, image_shift, box_geo);
    p.image_box() = p_ref.image_box() - image_shift;

    p.v() = velocity(p_ref, p.vs_relative());

    if (box_geo.type() == BoxType::LEES_EDWARDS) {
      auto const &shear_dir = box_geo.clees_edwards_bc().shear_direction;
      auto const &shear_normal = box_geo.clees_edwards_bc().shear_plane_normal;
      auto const &le_vel = box_geo.lees_edwards_bc().shear_velocity;
      Utils::Vector3i n_shifts{};
      fold_position(new_pos, n_shifts, box_geo);
      p.v()[shear_dir] -= n_shifts[shear_normal] * le_vel;
    }

    if (have_quaternions())
      p.quat() = p_ref.quat() * p.vs_relative().quat;
  }

  if (cell_structure.check_resort_required(particles, skin)) {
    cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  }
}

// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void VirtualSitesRelative::back_transfer_forces_and_torques() const {
  cell_structure.ghosts_reduce_forces();

  init_forces_ghosts(cell_structure.ghost_particles());

  // Iterate over all the particles in the local cells
  for (auto &p : cell_structure.local_particles()) {
    // We only care about virtual particles
    if (!p.is_virtual())
      continue;
    auto &p_ref = get_reference_particle(p.vs_relative());

    // Add forces and torques
    p_ref.force() += p.force();
    p_ref.torque() +=
        vector_product(connection_vector(p_ref, p.vs_relative()), p.force()) +
        p.torque();
  }
}

// Rigid body contribution to scalar pressure and pressure tensor
Utils::Matrix<double, 3, 3> VirtualSitesRelative::pressure_tensor() const {
  Utils::Matrix<double, 3, 3> pressure_tensor = {};

  for (auto &p : cell_structure.local_particles()) {
    if (!p.is_virtual())
      continue;

    auto const &p_ref = get_reference_particle(p.vs_relative());

    pressure_tensor += constraint_stress(p.force(), p_ref, p.vs_relative());
  }

  return pressure_tensor;
}
#endif // VIRTUAL_SITES_RELATIVE
