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
#include "system/System.hpp"
#include "cells.hpp"
#include "particle_node.hpp"

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
  mutable std::vector<double> contact_times;
  mutable std::vector<std::vector<bool>> contacts;
  mutable std::vector<std::vector<double>> first_contact_times;
  mutable std::vector<double> instantaneous_contact_times;

  explicit ContactTimes(std::vector<int> const &ids, std::vector<int> const &target_ids, double contact_threshold): 
  PidTimeObservable(ids,  target_ids, contact_threshold){
    if (this->ids().size() < 1) {throw std::runtime_error("At least 1 particle in ids is required");}
    if (this->target_ids.size() < 1) {throw std::runtime_error("At least 1 particle in target_ids is required");}
    if (contact_threshold < 0) {throw std::runtime_error("The contact threshold must be a positive number.");}
    auto const &system = System::get_system();
    // Initialize the contact time
    int n_ids = this->ids().size();
    int n_target_ids = target_ids.size();  
    double time = system.get_sim_time();
    this->first_contact_times = std::vector<std::vector<double>>(n_ids, std::vector<double>(n_target_ids, time));
    this->contacts = std::vector<std::vector<bool>>(n_ids, std::vector<bool>(n_target_ids, false));
    }
    
/**  Checks if `target` is an element in `vec`
 */
  bool is_target_in_vec(const std::vector<int>& vec, int target) const {
    return std::find(vec.begin(), vec.end(), target) != vec.end();
  }
  
/**  When particles with indexes `index1` and `index2` are not in contact, 
 *   update the contact map `contacts` and store a 0 contact time in the contact series `contact_times`
 */
  void update_contact_times_when_not_in_contact(double time, int index1, int index2) const{
    if (this->contacts[index1][index2]){ // index1 and index2 are now not in contact but they were before
      // # Calculate the total contact time
      auto first_contact_time = this->first_contact_times[index1][index2];
      auto contact_time = time - first_contact_time;
      this->contacts[index1][index2] = false;
      this->contact_times.push_back(contact_time);
    }
    else{this->contact_times.push_back(0);} // index1 and index2 were also not in contact before
  }

/**  Cleans up the series of contact times in memory
 */
  void clean_contact_times()const{
    this -> contact_times.clear();
    this -> instantaneous_contact_times.clear();}

/**  Returns the series of contact times stored in memory
 */
  std::vector<double> get_contact_times_series() const{return this->contact_times;}

/**  Returns the contact times for the current configuration
 */
  std::vector<double> get_instantaneous_contact_times() const{
    auto const &system = System::get_system();
    double time = system.get_sim_time();
    // Check for particles still in contact 
    auto ids1=this->ids();
    auto ids2=this->target_ids;
    for (auto &id1: ids1){
      for (auto &id2: ids2){
        auto index1=std::find(ids1.begin(), ids1.end(), id1)-ids1.begin();
        auto index2=std::find(ids2.begin(), ids2.end(), id2)-ids2.begin();
        if (this->contacts[index1][index2]){ // index1 and index2 are now in contact
          // # Calculate the current contact time
          auto first_contact_time = this->first_contact_times[index1][index2];
          auto contact_time = time - first_contact_time;
          this->instantaneous_contact_times.push_back(contact_time);
        }
      }
    }
  return this->instantaneous_contact_times;
  }

/**  Evaluates the current contact times
 */
  std::vector<double> evaluate(boost::mpi::communicator const &comm,
           ParticleReferenceRange const &local_particles,
           const ParticleObservables::traits<Particle> &traits) const override {

    if (comm.rank() != 0) {return {};}
    // Get instances of the system and cell structures
    auto const &system = System::get_system();
    auto const &box_geo  = *system.box_geo;
    auto &cell_structure = *system.cell_structure;
    auto neighbor_map = get_neighbor_pids(system);
    double time = system.get_sim_time();
    auto ids1=this->ids();
    auto ids2=this->target_ids;

    // Update the contact times
    for (auto &id1: ids1){
      auto neighbor_pids = neighbor_map[id1].neighbor_pids;
      for (auto &id2: ids2){
        // avoid counting contact times of a particle with itself
        if (id1 == id2){continue;}
        auto index1=std::find(ids1.begin(), ids1.end(), id1)-ids1.begin();
        auto index2=std::find(ids2.begin(), ids2.end(), id2)-ids2.begin();
        auto p1 = cell_structure.get_local_particle(id1);
        auto p2 = cell_structure.get_local_particle(id2);
        auto const dist =  box_geo.get_mi_vector(p1->pos(), p2->pos()).norm();  
        if (dist < contact_threshold) { // pid1 and pid2 are in contact now
          if (!(this->contacts[index1][index2])){ // but they were not in contact before!
              this->contacts[index1][index2]=true;
              this->first_contact_times[index1][index2]=time;
            }}
        else{update_contact_times_when_not_in_contact(time, index1,index2);} // // pid1 and pid2 are not in contact now
        
        
        
        /*
        if (is_target_in_vec(neighbor_pids,id2)){ // id1 and id2 are in the same neighbor list
          auto p1 = cell_structure.get_local_particle(id1);
          auto p2 = cell_structure.get_local_particle(id2);
          auto const dist =  box_geo.get_mi_vector(p1->pos(), p2->pos()).norm();  
          if (dist < contact_threshold) { // pid1 and pid2 are in contact now
            if (!(this->contacts[index1][index2])){ // but they were not in contact before!
              this->contacts[index1][index2]=true;
              this->first_contact_times[index1][index2]=time;
            }
          }
          else{update_contact_times_when_not_in_contact(time, index1,index2);} // // pid1 and pid2 are not in contact now               
        }
        else{update_contact_times_when_not_in_contact(time, index1,index2);} // // pid1 and pid2 are not in contact now
        */
      }
    }
    return {}; 
  }
  std::vector<std::size_t> shape_contact_time_series() const  {
    return {contact_times.size()};
  }
  std::vector<std::size_t> shape_instantaneous_contact_time() const {
    return {instantaneous_contact_times.size()};
  }
  std::vector<std::size_t> shape() const override {
    return {};
  }
};

} // namespace Observables

#endif
