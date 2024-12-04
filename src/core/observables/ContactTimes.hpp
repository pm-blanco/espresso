/*
 * Copyright (C) 2019-2022 The ESPResSo project
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
#ifndef OBSERVABLES_CONTACTTIMES_HPP
#define OBSERVABLES_CONTACTTIMES_HPP

#include "BoxGeometry.hpp"
#include "PidTimeObservable.hpp"
#include "system/System.hpp"

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>
#include <iostream>
#include "cell_system/CellStructure.hpp"

namespace Observables {

/** Calculates the contact times between the ids
 */
class ContactTimes : public PidTimeObservable {
public:
  using PidTimeObservable::PidTimeObservable;
  explicit ContactTimes(std::vector<int> const &ids, double max_z): 
PidTimeObservable(ids,  max_z){
  std::cout << "ContactTimes constructor called\n";  // Debugging line
    if (this->ids().size() < 2) {throw std::runtime_error("At least 2 particles are required");}
      
    if (max_z > 1) {std::cout<< max_z <<"\n";}
   } 

  std::vector<double>
  evaluate(boost::mpi::communicator const &comm,
           ParticleReferenceRange const &local_particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    auto const positions_sorted = detail::get_all_particle_positions(
        comm, local_particles, ids(), traits, false);
  
    if (comm.rank() != 0) {
      return {};
    }
    auto const &box_geo = *System::get_system().box_geo;
    std::vector<std::pair<Particle *, Particle *>> verlet_list = System::get_system().cell_structure->m_verlet_list;
    double time = System::get_system().get_sim_time();
    std::cout << time  << "\n";
    for (auto &particle_pair : verlet_list) {
      auto pos1 = particle_pair.first  -> pos();
      auto pos2 = particle_pair.second -> pos();
      auto const dist =  box_geo.get_mi_vector(pos1, pos2).norm();
      std::cout << max_z << "\n";
      if (dist < max_z) {
      std::cout << particle_pair.first -> id()  << " " << particle_pair.second -> id()  <<  " "  << dist << " " << max_z << "\n";
      }
    }
    
    std::vector<double> res(n_values());

    for (std::size_t i = 0, end = n_values(); i < end; i++) {
      auto const v =
          box_geo.get_mi_vector(positions_sorted[i], positions_sorted[i + 1]);
      res[i] = v.norm();
    }
    return res;
  }
  std::vector<std::size_t> shape() const override {
    assert(!ids().empty());
    return {ids().size() - 1};
  }
};

} // namespace Observables

#endif
