/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef OBSERVABLES_PARTICLE_TRAITS
#define OBSERVABLES_PARTICLE_TRAITS

#include "Particle.hpp"
#include "config.hpp"

namespace ParticleObservables {
/**
 * Template specialization for `Particle`. The traits mechanism is used to get
 * indirect access to particle properties. This helps making the implementation
 * of observables independent of the particle type.
 */
template <> struct traits<Particle> {
  auto position(Particle const &p) const { return p.pos(); }
  auto velocity(Particle const &p) const { return p.v(); }
  auto mass(Particle const &p) const {
#ifdef VIRTUAL_SITES
    // we exclude virtual particles since their mass does not have a meaning
    if (p.is_virtual())
      return decltype(p.mass()){};
#endif
    return p.mass();
  }
  auto charge(Particle const &p) const { return p.q(); }
  auto dipole_moment(Particle const &p) const {
#if defined(ROTATION) && defined(DIPOLES)
    return p.calc_dip();
#else
    return Utils::Vector3d{};
#endif
  }
};

} // namespace ParticleObservables

#endif
