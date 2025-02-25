/*
 * Copyright (C) 2025 The ESPResSo project
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
#include "cells.hpp"
#include "particle_node.hpp"
#include "grid.hpp"
#include "integrate.hpp"

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>
#include <iostream>
#include "cell_system/CellStructure.hpp"

namespace Observables {

/**  Tracks the time evolution of contacts between `ids` and `target_ids` within a given `contact_threshold`
 */
class ContactTimes : public PidTimeObservable {
public:
  using PidTimeObservable::PidTimeObservable;

  

  explicit ContactTimes(std::vector<int> const &ids, std::vector<int> const &target_ids, double contact_threshold): 
  PidTimeObservable(ids,  target_ids, contact_threshold),m_pairs{get_unique_pairs(ids, target_ids)}{
    if (this->ids().size() < 1) {throw std::runtime_error("At least 1 particle in ids is required");}
    if (this->target_ids.size() < 1) {throw std::runtime_error("At least 1 particle in target_ids is required");}
    if (contact_threshold < 0) {throw std::runtime_error("The contact threshold must be a positive number.");}
        
    }
    

/**  Cleans up the series of contact times in memory
 */
  void clean_contact_times()const{
    this -> contact_times.clear();}

/**  Returns the series of contact times stored in `contact_times`
 */
  std::vector<double> get_contact_times_series() const{return this->contact_times;}

/**  Evaluates the current contact times
 */
  std::vector<double> evaluate(ParticleReferenceRange particles,
           const ParticleObservables::traits<Particle> &) const override {

    double time = get_sim_time();
    double dt = get_time_step();
    // Initialize the bookkeeping vectors
    if (this->contacts.empty()) {
      this->contacts.resize(this->m_pairs.size(), false);
    }

    if (this->first_contact_times.empty()) {
      this->first_contact_times.resize(this->m_pairs.size(), time);
    }
    // Update the contact times
    auto index=0;
    for (const auto &p : m_pairs) {            
        auto p1 = cell_structure.get_local_particle(p.first);
        auto p2 = cell_structure.get_local_particle(p.second);
        auto const dist =  box_geo.get_mi_vector(p1->pos(), p2->pos()).norm();  
        
        if (dist < contact_threshold) { // pid1 and pid2 are in contact now
          if (!(this->contacts[index])){ // but they were not in contact before!
              this->contacts[index]=true;
              this->first_contact_times[index]=time;
            }}
        else{ // pid1 and pid2 are not in contact now
          if (this->contacts[index]){ // index1 and index2 are not in contact now but they were before
            // # Calculate the total contact time
            auto first_contact_time = this->first_contact_times[index];
            auto contact_time = time - first_contact_time - dt;
            if (contact_time < dt){contact_time=0;}
            this->contacts[index] = false;
            this->contact_times.push_back(contact_time);
          }
      
    } 
    index+=1; 
    }
    return {}; 
  }
  std::vector<std::size_t> shape_contact_time_series() const  {
    return {contact_times.size()};
  }
  
  std::vector<std::size_t> shape() const override {
    return {};
  }

  private:
    mutable std::vector<std::pair<int, int>> m_pairs;
    mutable std::vector<double> first_contact_times;
    mutable std::vector<double> contact_times;
    mutable std::vector<bool> contacts;
    mutable bool initialization;
    std::vector<std::pair<int, int>>
    get_unique_pairs(std::vector<int> const &ids1, std::vector<int> const &ids2) {
    std::set<std::pair<int, int>> unique_pairs;
    for (int id1 : ids1) {
      for (int id2 : ids2) {
        if (id1 != id2) {
          unique_pairs.emplace(std::minmax(id1, id2));
        }
      }
    }
    return {unique_pairs.begin(), unique_pairs.end()};
  }

};

} // namespace Observables

#endif
