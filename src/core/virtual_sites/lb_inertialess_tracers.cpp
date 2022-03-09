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
/// \file
/// \brief Main of the Bayreuth Immersed-Boundary implementation

#include "config.hpp"

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS

#include "virtual_sites/lb_inertialess_tracers.hpp"

#include "Particle.hpp"
#include "cells.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "integrate.hpp"
#include "lb_inertialess_tracers_cuda_interface.hpp"

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

#include <cstddef>

void CoupleIBMParticleToFluid(Particle const &p);
void ParticleVelocitiesFromLB_CPU();
bool IsHalo(std::size_t indexCheck);

static bool *isHaloCache = nullptr;

inline bool in_local_domain(Utils::Vector3d const &pos) {
  auto const offset = Utils::Vector3d::broadcast(0.5 * lblattice.agrid);
  auto const my_left = local_geo.my_left();
  auto const my_right = local_geo.my_right();

  return pos >= (my_left - offset) and pos < (my_right + offset);
}

/** Put the calculated force stored on the ibm particles into the fluid by
 *  updating the @ref lbfields structure.
 *  Called from the integration loop right after the forces have been
 *  calculated.
 */
void IBM_ForcesIntoFluid_CPU() {
  // Update the forces on the ghost particles
  cell_structure.ghosts_update(Cells::DATA_PART_FORCE);

  // Loop over local cells
  for (auto &p : cell_structure.local_particles()) {
    if (p.is_virtual())
      CoupleIBMParticleToFluid(p);
  }

  for (auto &p : cell_structure.ghost_particles()) {
    // for ghost particles we have to check if they lie
    // in the range of the local lattice nodes
    if (in_local_domain(p.pos())) {
      if (p.is_virtual())
        CoupleIBMParticleToFluid(p);
    }
  }
}

/** Interpolate LB velocity at the particle positions and propagate the
 *  particles.
 *  Called from the integration loop right after the LB update.
 */
void IBM_UpdateParticlePositions(ParticleRange const &particles,
                                 double time_step, int this_node) {
  // Get velocities
  if (lattice_switch == ActiveLB::CPU)
    ParticleVelocitiesFromLB_CPU();
#ifdef CUDA
  if (lattice_switch == ActiveLB::GPU)
    ParticleVelocitiesFromLB_GPU(particles, this_node);
#endif

  // Euler integrator
  for (auto &p : particles) {
    if (p.is_virtual()) {
      for (int axis = 0; axis < 2; axis++) {
#ifdef EXTERNAL_FORCES
        if (not p.is_fixed_along(axis))
#endif
          p.pos()[axis] += p.v()[axis] * time_step;
      }
    }
  }

  if (cell_structure.check_resort_required(particles, skin)) {
    cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  }
}

/** Put the momentum of a given particle into the LB fluid. */
void CoupleIBMParticleToFluid(Particle const &p) {
  // Convert units from MD to LB
  auto const delta_j = p.force() * Utils::int_pow<4>(lbpar.tau) / lbpar.agrid;

  // Get indices and weights of affected nodes using discrete delta function
  Utils::Vector<std::size_t, 8> node_index{};
  Utils::Vector6d delta{};
  lblattice.map_position_to_lattice(p.pos(), node_index, delta);

  // Loop over all affected nodes
  for (int z = 0; z < 2; z++) {
    for (int y = 0; y < 2; y++) {
      for (int x = 0; x < 2; x++) {
        // Do not put force into a halo node
        if (!IsHalo(static_cast<int>(node_index[(z * 2 + y) * 2 + x]))) {
          // Add force into the lbfields structure
          auto &local_f =
              lbfields[node_index[(z * 2 + y) * 2 + x]].force_density;

          local_f +=
              delta[3 * x + 0] * delta[3 * y + 1] * delta[3 * z + 2] * delta_j;
        }
      }
    }
  }
}

/** Calculate the LB fluid velocity at a particle position.
 *  Very similar to the velocity interpolation done in standard ESPResSo,
 *  except that we add the f/2 contribution, cf. @cite guo02a.
 *  The fluid velocity is obtained by linear interpolation,
 *  cf. eq. (11) in @cite ahlrichs99a.
 */
template <bool ReturnVelocity>
Utils::Vector3d GetIBMInterpolatedVelocity(Utils::Vector3d const &pos) {
  auto const f_ext =
      lbpar.ext_force_density * Utils::sqr(lbpar.agrid * lbpar.tau);

  /* determine elementary lattice cell surrounding the particle
     and the relative position of the particle in this cell */
  Utils::Vector<std::size_t, 8> node_index{};
  Utils::Vector6d delta{};
  lblattice.map_position_to_lattice(pos, node_index, delta);

  // This for the f/2 contribution to the velocity
  Utils::Vector3d force_added = {};
  Utils::Vector3d interpolated_u = {};

  for (int z = 0; z < 2; z++) {
    for (int y = 0; y < 2; y++) {
      for (int x = 0; x < 2; x++) {
        auto const index = node_index[(z * 2 + y) * 2 + x];
        auto const local_delta =
            delta[3 * x + 0] * delta[3 * y + 1] * delta[3 * z + 2];
        const auto &f = lbfields[index].force_density_buf;

        double local_density;
        Utils::Vector3d local_j;

        // This can be done more easily without copying the code twice.
        // We probably can even set the boundary velocity directly.
#ifdef LB_BOUNDARIES
        if (lbfields[index].boundary) {
          if (ReturnVelocity) {
            local_density = lbpar.density;
            auto const i = lbfields[index].boundary - 1;
            local_j = lbpar.density * LBBoundaries::lbboundaries[i]->velocity();
          }
        } else
#endif
        {
          auto const modes = lb_calc_modes(static_cast<int>(index), lbfluid);
          local_density = lbpar.density + modes[0];

          if (ReturnVelocity) {
            // Add the +f/2 contribution!!
            local_j[0] = modes[1] + f[0] / 2.;
            local_j[1] = modes[2] + f[1] / 2.;
            local_j[2] = modes[3] + f[2] / 2.;
          } else {
            // Keep track of the forces that we added to the fluid.
            // This is necessary for communication because this part is executed
            // for real and ghost particles.
            // Later on we sum the real and ghost contributions.
            force_added += local_delta * (f - f_ext) / (2. * local_density);
          }
        }

        // Interpolate velocity
        if (ReturnVelocity) {
          interpolated_u += local_j * (local_delta / local_density);
        }
      }
    }
  }

  auto const unit_conversion = lbpar.agrid / lbpar.tau;
  if (ReturnVelocity) {
    return interpolated_u * unit_conversion;
  }
  return force_added * unit_conversion;
}

/** Build a cache structure which contains a flag for each LB node whether
 * that node is a halo node or not.
 */
bool IsHalo(std::size_t indexCheck) {
  // First call --> build cache
  if (isHaloCache == nullptr) {
    isHaloCache = new bool[lblattice.halo_grid_volume];
    // Assume everything is a halo and correct in the next step
    for (int i = 0; i < lblattice.halo_grid_volume; i++)
      isHaloCache[i] = true;
    // Loop through and check where indexCheck occurs
    auto index = lblattice.halo_offset;
    for (int z = 1; z <= lblattice.grid[2]; z++) {
      for (int y = 1; y <= lblattice.grid[1]; y++) {
        for (int x = 1; x <= lblattice.grid[0]; x++) {
          isHaloCache[index] = false;
          ++index;
        }
        index += 2; /* skip halo region */
      }
      index += 2 * lblattice.halo_grid[0]; /* skip halo region */
    }
  }

  // Return
  return isHaloCache[indexCheck];
}

/** Get particle velocities from LB and set the velocity field in the
 * particles data structure.
 */
void ParticleVelocitiesFromLB_CPU() {
  // Loop over particles in local cells.
  // Here all contributions are included: velocity, external force and
  // particle force.
  for (auto &p : cell_structure.local_particles()) {
    if (p.is_virtual()) {
      // Get interpolated velocity and store in the force (!) field
      // for later communication (see below)
      p.force() = GetIBMInterpolatedVelocity<true>(p.pos());
    }
  }

  // Loop over particles in ghost cells
  // Here we only add the particle forces stemming from the ghosts
  for (auto &p : cell_structure.ghost_particles()) {
    // This criterion include the halo on the left, but excludes the halo on
    // the right
    // Try if we have to use *1.5 on the right
    if (in_local_domain(p.pos())) {
      if (p.is_virtual()) {
        // The force stemming from the ghost particle (for
        // communication, see below)
        p.force() = GetIBMInterpolatedVelocity<false>(p.pos());
      } else {
        // Reset, necessary because we add all forces below. Also needs to
        // be done for real particles!
        p.force() = {};
      }
    } else {
      // Reset, necessary because we add all forces below
      p.force() = {};
    }
  }

  // Now the local particles contain a velocity (stored in the force field)
  // and the ghosts contain the rest of the velocity in their respective force
  // fields.
  // We need to add these. Since we have stored them in the force, not the
  // velocity fields, we can use the standard force communicator and then
  // transfer to the velocity afterwards.
  // Note that this overwrites the actual force which would be a problem for
  // real particles.
  // This could be solved by keeping a backup of the local forces before this
  // operation is attempted.
  cell_structure.ghosts_reduce_forces();

  // Transfer to velocity field
  for (auto &p : cell_structure.local_particles()) {
    if (p.is_virtual()) {
      p.v() = p.force();
    }
  }
}
#endif
