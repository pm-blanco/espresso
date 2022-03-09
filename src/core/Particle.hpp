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
#ifndef ESPRESSO_CORE_PARTICLE_HPP
#define ESPRESSO_CORE_PARTICLE_HPP

#include "config.hpp"

#include "BondList.hpp"

#include <utils/Vector.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/quaternion.hpp>

#include <boost/serialization/vector.hpp>

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace detail {
inline void check_axis_idx_valid(int const axis) {
  assert(axis >= 0 and axis <= 2);
}

inline bool get_nth_bit(uint8_t const bitfield, int const bit_idx) {
  return bitfield & (1u << bit_idx);
}
} // namespace detail

/** Properties of a self-propelled particle. */
struct ParticleParametersSwimming {
  /** Is the particle a swimmer. */
  bool swimming = false;
  /** Imposed constant force. */
  double f_swim = 0.;
  /** Constant velocity to relax to. */
  double v_swim = 0.;
  /** Flag for the swimming mode in a LB fluid.
   *  Values:
   *  - -1: pusher
   *  - +1: puller
   *  - 0: no swimming
   */
  int push_pull = 0;
  /** Distance of the source of propulsion from the particle
   *  center in a LB fluid.
   */
  double dipole_length = 0.;

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &swimming &f_swim &v_swim &push_pull &dipole_length;
  }
};

/** Properties of a particle which are not supposed to
 *  change during the integration, but have to be known
 *  for all ghosts. Ghosts are particles which are
 *  needed in the interaction calculation, but are just copies of
 *  particles stored on different nodes.
 */
struct ParticleProperties {
  /** unique identifier for the particle. */
  int identity = -1;
  /** Molecule identifier. */
  int mol_id = 0;
  /** particle type, used for non-bonded interactions. */
  int type = 0;

  /** particle mass */
#ifdef MASS
  double mass = 1.0;
#else
  constexpr static double mass{1.0};
#endif

  /** rotational inertia */
#ifdef ROTATIONAL_INERTIA
  Utils::Vector3d rinertia = {1., 1., 1.};
#else
  static constexpr Utils::Vector3d rinertia = {1., 1., 1.};
#endif

#ifdef ROTATION
  /** Bitfield for the particle axes of rotation.
   *  Values:
   *  - 1: allow rotation around the x axis
   *  - 2: allow rotation around the y axis
   *  - 4: allow rotation around the z axis
   *  By default, the particle cannot rotate.
   */
  uint8_t rotation = static_cast<uint8_t>(0b000u);
#else
  /** Bitfield for the particle axes of rotation. Particle cannot rotate. */
  static constexpr uint8_t rotation = static_cast<uint8_t>(0b000u);
#endif

  /** charge. */
#ifdef ELECTROSTATICS
  double q = 0.0;
#else
  constexpr static double q{0.0};
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
  /** electrophoretic mobility times E-field: mu_0 * E */
  Utils::Vector3d mu_E = {0., 0., 0.};
#endif

#ifdef DIPOLES
  /** dipole moment (absolute value) */
  double dipm = 0.;
#endif

#ifdef VIRTUAL_SITES
  /** is particle virtual */
  bool is_virtual = false;
#ifdef VIRTUAL_SITES_RELATIVE
  /** The following properties define, with respect to which real particle a
   *  virtual site is placed and at what distance. The relative orientation of
   *  the vector pointing from real particle to virtual site with respect to the
   *  orientation of the real particle is stored in the virtual site's
   *  quaternion attribute.
   */
  struct VirtualSitesRelativeParameters {
    int to_particle_id = 0;
    double distance = 0;
    /** Relative position of the virtual site. */
    Utils::Quaternion<double> rel_orientation =
        Utils::Quaternion<double>::identity();
    /** Orientation of the virtual particle in the body fixed frame. */
    Utils::Quaternion<double> quat = Utils::Quaternion<double>::identity();

    template <class Archive> void serialize(Archive &ar, long int) {
      ar &to_particle_id;
      ar &distance;
      ar &rel_orientation;
      ar &quat;
    }
  } vs_relative;
#endif // VIRTUAL_SITES_RELATIVE
#else  // VIRTUAL_SITES
  static constexpr bool is_virtual = false;
#endif // VIRTUAL_SITES

#ifdef THERMOSTAT_PER_PARTICLE
/** Friction coefficient for translation */
#ifndef PARTICLE_ANISOTROPY
  double gamma = -1.;
#else
  Utils::Vector3d gamma = {-1., -1., -1.};
#endif // PARTICLE_ANISOTROPY
#ifdef ROTATION
/** Friction coefficient for rotation */
#ifndef PARTICLE_ANISOTROPY
  double gamma_rot = -1.;
#else
  Utils::Vector3d gamma_rot = {-1., -1., -1.};
#endif // PARTICLE_ANISOTROPY
#endif // ROTATION
#endif // THERMOSTAT_PER_PARTICLE

#ifdef EXTERNAL_FORCES
  /** Bitfield for fixed particle coordinates.
   *  Values:
   *  - 1: fix translation along the x axis
   *  - 2: fix translation along the y axis
   *  - 4: fix translation along the z axis
   */
  uint8_t ext_flag = static_cast<uint8_t>(0b000u);
  /** External force. */
  Utils::Vector3d ext_force = {0, 0, 0};
#ifdef ROTATION
  /** External torque. */
  Utils::Vector3d ext_torque = {0, 0, 0};
#endif
#else  // EXTERNAL_FORCES
  /** Bitfield for fixed particle coordinates. Coordinates cannot be fixed. */
  static constexpr const uint8_t ext_flag = static_cast<uint8_t>(0b000u);
#endif // EXTERNAL_FORCES

#ifdef ENGINE
  ParticleParametersSwimming swim;
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &identity;
    ar &mol_id;
    ar &type;
#ifdef MASS
    ar &mass;
#endif
#ifdef ROTATIONAL_INERTIA
    ar &rinertia;
#endif
#ifdef ROTATION
    ar &rotation;
#endif
#ifdef ELECTROSTATICS
    ar &q;
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
    ar &mu_E;
#endif
#ifdef DIPOLES
    ar &dipm;
#endif

#ifdef VIRTUAL_SITES
    ar &is_virtual;
#ifdef VIRTUAL_SITES_RELATIVE
    ar &vs_relative;
#endif
#endif // VIRTUAL_SITES

#ifdef THERMOSTAT_PER_PARTICLE
    ar &gamma;
#ifdef ROTATION
    ar &gamma_rot;
#endif
#endif // THERMOSTAT_PER_PARTICLE
#ifdef EXTERNAL_FORCES
    ar &ext_flag;
    ar &ext_force;
#ifdef ROTATION
    ar &ext_torque;
#endif
#endif // EXTERNAL_FORCES

#ifdef ENGINE
    ar &swim;
#endif
  }
};

/** Positional information on a particle. Information that is
 *  communicated to calculate interactions with ghost particles.
 */
struct ParticlePosition {
  /** periodically folded position. */
  Utils::Vector3d p = {0, 0, 0};

#ifdef ROTATION
  /** quaternion to define particle orientation */
  Utils::Quaternion<double> quat = Utils::Quaternion<double>::identity();
  /** unit director calculated from the quaternion */
  Utils::Vector3d calc_director() const {
    return Utils::convert_quaternion_to_director(quat);
  };
#endif

#ifdef BOND_CONSTRAINT
  /** particle position at the previous time step (RATTLE algorithm) */
  Utils::Vector3d p_last_timestep = {0., 0., 0.};
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &p;
#ifdef ROTATION
    ar &quat;
#endif
#ifdef BOND_CONSTRAINT
    ar &p_last_timestep;
#endif
  }
};

/** Force information on a particle. Forces of ghost particles are
 *  collected and added up to the force of the original particle.
 */
struct ParticleForce {
  ParticleForce() = default;
  ParticleForce(ParticleForce const &) = default;
  ParticleForce &operator=(ParticleForce const &) = default;
  ParticleForce(const Utils::Vector3d &f) : f(f) {}
#ifdef ROTATION
  ParticleForce(const Utils::Vector3d &f, const Utils::Vector3d &torque)
      : f(f), torque(torque) {}
#endif

  friend ParticleForce operator+(ParticleForce const &lhs,
                                 ParticleForce const &rhs) {
#ifdef ROTATION
    return {lhs.f + rhs.f, lhs.torque + rhs.torque};
#else
    return lhs.f + rhs.f;
#endif
  }

  ParticleForce &operator+=(ParticleForce const &rhs) {
    return *this = *this + rhs;
  }

  /** force. */
  Utils::Vector3d f = {0., 0., 0.};

#ifdef ROTATION
  /** torque. */
  Utils::Vector3d torque = {0., 0., 0.};
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &f;
#ifdef ROTATION
    ar &torque;
#endif
  }
};

/** Momentum information on a particle. Information not contained in
 *  communication of ghost particles so far, but a communication would
 *  be necessary for velocity-dependent potentials.
 */
struct ParticleMomentum {
  /** velocity. */
  Utils::Vector3d v = {0., 0., 0.};

#ifdef ROTATION
  /** angular velocity.
   *  ALWAYS IN PARTICLE FIXED, I.E., CO-ROTATING COORDINATE SYSTEM.
   */
  Utils::Vector3d omega = {0., 0., 0.};
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &v;
#ifdef ROTATION
    ar &omega;
#endif
  }
};

/** Information on a particle that is needed only on the
 *  node the particle belongs to.
 */
struct ParticleLocal {
  /** is particle a ghost particle. */
  bool ghost = false;
  /** position from the last Verlet list update. */
  Utils::Vector3d p_old = {0, 0, 0};
  /** index of the simulation box image where the particle really sits. */
  Utils::Vector3i i = {0, 0, 0};

  /** Accumulated applied Lees-Edwards offset. */
  double lees_edwards_offset = 0.;
  short int lees_edwards_flag = 0;

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &ghost;
    ar &p_old;
    ar &i;
    ar &lees_edwards_offset;
    ar &lees_edwards_flag;
  }
};

#ifdef BOND_CONSTRAINT
struct ParticleRattle {
  /** position/velocity correction */
  Utils::Vector3d correction = {0, 0, 0};

  friend ParticleRattle operator+(ParticleRattle const &lhs,
                                  ParticleRattle const &rhs) {
    return {lhs.correction + rhs.correction};
  }

  ParticleRattle &operator+=(ParticleRattle const &rhs) {
    return *this = *this + rhs;
  }

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &correction;
  }
};
#endif

/** Struct holding all information for one particle. */
struct Particle { // NOLINT(bugprone-exception-escape)
  int &identity() { return p.identity; }
  int const &identity() const { return p.identity; }

  bool operator==(Particle const &rhs) const {
    return identity() == rhs.identity();
  }

  bool operator!=(Particle const &rhs) const {
    return identity() != rhs.identity();
  }

  ///
  ParticleProperties p;
  ///
  ParticlePosition r;
#ifdef DIPOLES
  Utils::Vector3d calc_dip() const { return r.calc_director() * p.dipm; }
#endif
  ///
  ParticleMomentum m;
  ///
  ParticleForce f;
  ///
  ParticleLocal l;
#ifdef BOND_CONSTRAINT
  ///
  ParticleRattle rattle;
#endif

private:
  BondList bl;

#ifdef EXCLUSIONS
  /** list of particles, with which this particle has no non-bonded
   *  interactions
   */
  std::vector<int> el;
#endif

public:
  auto const &id() const { return p.identity; }
  auto &id() { return p.identity; }
  auto const &mol_id() const { return p.mol_id; }
  auto &mol_id() { return p.mol_id; }
  auto const &type() const { return p.type; }
  auto &type() { return p.type; }

  auto const &bonds() const { return bl; }
  auto &bonds() { return bl; }

  auto const &pos() const { return r.p; }
  auto &pos() { return r.p; }
  auto const &v() const { return m.v; }
  auto &v() { return m.v; }
  auto const &force() const { return f.f; }
  auto &force() { return f.f; }

  bool is_ghost() const { return l.ghost; }
  void set_ghost(bool const ghost_flag) { l.ghost = ghost_flag; }
  auto &pos_at_last_verlet_update() { return l.p_old; }
  auto const &pos_at_last_verlet_update() const { return l.p_old; }
  auto const &image_box() const { return l.i; }
  auto &image_box() { return l.i; }
  auto const &lees_edwards_offset() const { return l.lees_edwards_offset; }
  auto &lees_edwards_offset() { return l.lees_edwards_offset; }
  auto const &lees_edwards_flag() const { return l.lees_edwards_flag; }
  auto &lees_edwards_flag() { return l.lees_edwards_flag; }

#ifdef MASS
  auto const &mass() const { return p.mass; }
  auto &mass() { return p.mass; }
#else
  constexpr auto &mass() const { return p.mass; }
#endif
#ifdef ROTATION
  bool can_rotate() const { return static_cast<bool>(p.rotation); }
  bool can_rotate_around(int const axis) const {
    detail::check_axis_idx_valid(axis);
    return detail::get_nth_bit(p.rotation, axis);
  }
  void set_can_rotate_around(int const axis, bool const rot_flag) {
    detail::check_axis_idx_valid(axis);
    if (rot_flag) {
      p.rotation |= static_cast<uint8_t>(1u << axis);
    } else {
      p.rotation &= static_cast<uint8_t>(~(1u << axis));
    }
  }
  void set_can_rotate_all_axes() { p.rotation = static_cast<uint8_t>(0b111u); }
  void set_cannot_rotate_all_axes() {
    p.rotation = static_cast<uint8_t>(0b000u);
  }
  auto const &quat() const { return r.quat; }
  auto &quat() { return r.quat; }
  auto const &torque() const { return f.torque; }
  auto &torque() { return f.torque; }
  auto const &omega() const { return m.omega; }
  auto &omega() { return m.omega; }
  auto const &ext_torque() const { return p.ext_torque; }
  auto &ext_torque() { return p.ext_torque; }
  auto calc_director() const { return r.calc_director(); }
#else
  bool can_rotate() const { return false; }
  bool can_rotate_around(int const axis) const { return false; }
#endif
#ifdef DIPOLES
  auto const &dipm() const { return p.dipm; }
  auto &dipm() { return p.dipm; }
#endif
#ifdef ROTATIONAL_INERTIA
  auto const &rinertia() const { return p.rinertia; }
  auto &rinertia() { return p.rinertia; }
#else
  constexpr auto &rinertia() const { return p.rinertia; }
#endif
#ifdef ELECTROSTATICS
  auto const &q() const { return p.q; }
  auto &q() { return p.q; }
#else
  constexpr auto &q() const { return p.q; }
#endif
#ifdef LB_ELECTROHYDRODYNAMICS
  auto const &mu_E() const { return p.mu_E; }
  auto &mu_E() { return p.mu_E; }
#endif
#ifdef VIRTUAL_SITES
  bool is_virtual() const { return p.is_virtual; }
  void set_virtual(bool const virt_flag) { p.is_virtual = virt_flag; }
#ifdef VIRTUAL_SITES_RELATIVE
  auto const &vs_relative() const { return p.vs_relative; }
  auto &vs_relative() { return p.vs_relative; }
#endif // VIRTUAL_SITES_RELATIVE
#else
  constexpr auto is_virtual() const { return p.is_virtual; }
#endif
#ifdef THERMOSTAT_PER_PARTICLE
  auto const &gamma() const { return p.gamma; }
  auto &gamma() { return p.gamma; }
#ifdef ROTATION
  auto const &gamma_rot() const { return p.gamma_rot; }
  auto &gamma_rot() { return p.gamma_rot; }
#endif // ROTATION
#endif // THERMOSTAT_PER_PARTICLE
#ifdef EXTERNAL_FORCES
  bool has_fixed_coordinates() const { return static_cast<bool>(p.ext_flag); }
  bool is_fixed_along(int const axis) const {
    detail::check_axis_idx_valid(axis);
    return detail::get_nth_bit(p.ext_flag, axis);
  }
  void set_fixed_along(int const axis, bool const fixed_flag) {
    // set new flag
    if (fixed_flag) {
      p.ext_flag |= static_cast<uint8_t>(1u << axis);
    } else {
      p.ext_flag &= static_cast<uint8_t>(~(1u << axis));
    }
  }
  auto const &ext_force() const { return p.ext_force; }
  auto &ext_force() { return p.ext_force; }

#else  // EXTERNAL_FORCES
  constexpr bool has_fixed_coordinates() const { return false; }
  constexpr bool is_fixed_along(int const) const { return false; }
#endif // EXTERNAL_FORCES
#ifdef ENGINE
  auto const &swimming() const { return p.swim; }
  auto &swimming() { return p.swim; }
#endif
#ifdef BOND_CONSTRAINT
  auto const &pos_last_time_step() const { return r.p_last_timestep; }
  auto &pos_last_time_step() { return r.p_last_timestep; }
  auto const &rattle_params() const { return rattle; }
  auto &rattle_params() { return rattle; }
#endif
#ifdef EXCLUSIONS
  auto const &exclusions() const { return el; }
  auto &exclusions() { return el; }
#endif

private:
  friend boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &p;
    ar &r;
    ar &m;
    ar &f;
    ar &l;
    ar &bl;
#ifdef EXCLUSIONS
    ar &el;
#endif
  }
};

#endif
