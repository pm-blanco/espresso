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

#include "reaction_ensemble.hpp"
#include "Particle.hpp"
#include "config.hpp"
#include "energy.hpp"
#include "grid.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "statistics.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/contains.hpp>
#include <utils/index.hpp>

#include <boost/range/algorithm.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ReactionEnsemble {

/** Load minimum and maximum energies as a function of the other collective
 *  variables.
 */
void EnergyCollectiveVariable::load_CV_boundaries(
    WangLandauReactionEnsemble &wl_system, std::istream &infile) {

  wl_system.do_energy_reweighting = true;

  // Note that you cannot change the other collective variables in the
  // pre-production run and the production run

  // The data is formatted as space-separated floaing point values
  // (the first line is a header). The min and max energies are stored in
  // the last two columns. The first N columns are the collective variables.
  std::string line = "";
  std::getline(infile, line); // dummy read to throw away header
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::vector<double> values;
    double value = -1.0;
    while (iss >> value) {
      values.push_back(value);
    }
    assert(values.size() >= 2);
    wl_system.max_boundaries_energies.emplace_back(values.back());
    wl_system.min_boundaries_energies.emplace_back(values[values.size() - 2]);
  }

  CV_minimum = *boost::range::min_element(wl_system.min_boundaries_energies);
  CV_maximum = *boost::range::max_element(wl_system.max_boundaries_energies);
}

/**
 * Performs a randomly selected reaction in the reaction ensemble
 */
int ReactionAlgorithm::do_reaction(int reaction_steps) {
  for (int i = 0; i < reaction_steps; i++) {
    int reaction_id = i_random(static_cast<int>(reactions.size()));
    generic_oneway_reaction(reaction_id);
  }
  return 0;
}

/**
 * Adds a reaction to the reaction system
 */
void ReactionAlgorithm::add_reaction(
    double gamma, const std::vector<int> &reactant_types,
    const std::vector<int> &reactant_coefficients,
    const std::vector<int> &product_types,
    const std::vector<int> &product_coefficients) {
  SingleReaction new_reaction(gamma, reactant_types, reactant_coefficients,
                              product_types, product_coefficients);

  // make ESPResSo count the particle numbers which take part in the reactions
  for (int reactant_type : new_reaction.reactant_types)
    init_type_map(reactant_type);
  for (int product_type : new_reaction.product_types)
    init_type_map(product_type);

  init_type_map(non_interacting_type);

  reactions.push_back(new_reaction);
}

/**
 * Checks whether all necessary variables for the reaction ensemble have been
 * set.
 */
void ReactionAlgorithm::check_reaction_ensemble() const {
  if (reactions.empty()) {
    throw std::runtime_error("Reaction system not initialized");
  }

  if (temperature < 0) {
    throw std::runtime_error("Temperatures cannot be negative. Please provide "
                             "a temperature (in k_B T) to the simulation. "
                             "Normally it should be 1.0. This will be used "
                             "directly to calculate beta:=1/(k_B T) which "
                             "occurs in the exp(-beta*E)\n");
  }
#ifdef ELECTROSTATICS
  // check for the existence of default charges for all types that take part in
  // the reactions

  for (const auto &current_reaction : reactions) {
    // check for reactants
    for (int reactant_type : current_reaction.reactant_types) {
      auto it = charges_of_types.find(reactant_type);
      if (it == charges_of_types.end()) {
        std::string message = std::string("Forgot to assign charge to type ") +
                              std::to_string(reactant_type);
        throw std::runtime_error(message);
      }
    }
    // check for products
    for (int product_type : current_reaction.product_types) {
      auto it = charges_of_types.find(product_type);
      if (it == charges_of_types.end()) {
        std::string message = std::string("Forgot to assign charge to type ") +
                              std::to_string(product_type);
        throw std::runtime_error(message);
      }
    }
  }
#endif
}

/**
 * Automatically sets the volume which is used by the reaction ensemble to the
 * volume of a cuboid box
 */
void ReactionAlgorithm::set_cuboid_reaction_ensemble_volume() {
  if (volume < 0)
    volume = box_geo.volume();
}

/**
 * Calculates the factorial expression which occurs in the reaction ensemble
 * acceptance probability
 */
double factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(int Ni0, int nu_i) {
  double value = 1.0;
  if (nu_i) {
    if (nu_i > 0) {
      for (int i = 1; i <= nu_i; i++) {
        value /= Ni0 + i;
      }
    } else {
      auto const abs_nu_i = std::abs(nu_i);
      for (int i = 0; i < abs_nu_i; i++) {
        value *= Ni0 - i;
      }
    }
  }
  return value;
}

/**
 * Checks whether all particles exist for the provided reaction.
 */
bool ReactionAlgorithm::all_reactant_particles_exist(int reaction_id) const {
  bool enough_particles = true;
  for (int i = 0; i < reactions[reaction_id].reactant_types.size(); i++) {
    int current_number =
        number_of_particles_with_type(reactions[reaction_id].reactant_types[i]);
    if (current_number < reactions[reaction_id].reactant_coefficients[i]) {
      enough_particles = false;
      break;
    }
  }
  return enough_particles;
}

/**
 * Stores the particle property of a random particle of the provided type into
 * the provided vector
 */
void ReactionAlgorithm::append_particle_property_of_random_particle(
    int type, std::vector<StoredParticleProperty> &list_of_particles) {
  int random_index_in_type_map = i_random(number_of_particles_with_type(type));
  int p_id = get_random_p_id(type, random_index_in_type_map);
  StoredParticleProperty property_of_part = {p_id, charges_of_types[type],
                                             type};
  list_of_particles.push_back(property_of_part);
}

/**
 *Performs a trial reaction move
 */
void ReactionAlgorithm::make_reaction_attempt(
    SingleReaction const &current_reaction,
    std::vector<StoredParticleProperty> &changed_particles_properties,
    std::vector<int> &p_ids_created_particles,
    std::vector<StoredParticleProperty> &hidden_particles_properties) {
  // create or hide particles of types with corresponding types in reaction
  for (int i = 0; i < std::min(current_reaction.product_types.size(),
                               current_reaction.reactant_types.size());
       i++) {
    // change std::min(reactant_coefficients(i),product_coefficients(i)) many
    // particles of reactant_types(i) to product_types(i)
    for (int j = 0; j < std::min(current_reaction.product_coefficients[i],
                                 current_reaction.reactant_coefficients[i]);
         j++) {
      append_particle_property_of_random_particle(
          current_reaction.reactant_types[i], changed_particles_properties);
      replace_particle(changed_particles_properties.back().p_id,
                       current_reaction.product_types[i]);
    }
    // create product_coefficients(i)-reactant_coefficients(i) many product
    // particles iff product_coefficients(i)-reactant_coefficients(i)>0,
    // iff product_coefficients(i)-reactant_coefficients(i)<0, hide this number
    // of reactant particles
    if (current_reaction.product_coefficients[i] -
            current_reaction.reactant_coefficients[i] >
        0) {
      for (int j = 0; j < current_reaction.product_coefficients[i] -
                              current_reaction.reactant_coefficients[i];
           j++) {
        int p_id = create_particle(current_reaction.product_types[i]);
        p_ids_created_particles.push_back(p_id);
      }
    } else if (current_reaction.reactant_coefficients[i] -
                   current_reaction.product_coefficients[i] >
               0) {
      for (int j = 0; j < current_reaction.reactant_coefficients[i] -
                              current_reaction.product_coefficients[i];
           j++) {
        append_particle_property_of_random_particle(
            current_reaction.reactant_types[i], hidden_particles_properties);
        hide_particle(hidden_particles_properties.back().p_id,
                      current_reaction.reactant_types[i]);
      }
    }
  }
  // create or hide particles of types with noncorresponding replacement types
  for (auto i = std::min(current_reaction.product_types.size(),
                         current_reaction.reactant_types.size());
       i < std::max(current_reaction.product_types.size(),
                    current_reaction.reactant_types.size());
       i++) {
    if (current_reaction.product_types.size() <
        current_reaction.reactant_types.size()) {
      // hide superfluous reactant_types particles
      for (int j = 0; j < current_reaction.reactant_coefficients[i]; j++) {
        append_particle_property_of_random_particle(
            current_reaction.reactant_types[i], hidden_particles_properties);
        hide_particle(hidden_particles_properties.back().p_id,
                      current_reaction.reactant_types[i]);
      }
    } else {
      // create additional product_types particles
      for (int j = 0; j < current_reaction.product_coefficients[i]; j++) {
        int p_id = create_particle(current_reaction.product_types[i]);
        p_ids_created_particles.push_back(p_id);
      }
    }
  }
}

/**
 * Restores the previously stored particle properties. This function is invoked
 * when a reaction attempt is rejected.
 */
void ReactionAlgorithm::restore_properties(
    std::vector<StoredParticleProperty> const &property_list,
    const int number_of_saved_properties) {
  // this function restores all properties of all particles provided in the
  // property list, the format of the property list is (p_id,charge,type)
  // repeated for each particle that occurs in that list
  for (auto &i : property_list) {
    int type = i.type;
#ifdef ELECTROSTATICS
    // set charge
    double charge = i.charge;
    set_particle_q(i.p_id, charge);
#endif
    // set type
    set_particle_type(i.p_id, type);
  }
}

/**
 * Calculates the expression in the acceptance probability in the reaction
 * ensemble
 */
double ReactionEnsemble::calculate_acceptance_probability(
    SingleReaction const &current_reaction, double E_pot_old, double E_pot_new,
    std::map<int, int> const &old_particle_numbers, int dummy_old_state_index,
    int dummy_new_state_index,
    bool dummy_only_make_configuration_changing_move) const {
  const double factorial_expr =
      calculate_factorial_expression(current_reaction, old_particle_numbers);

  const double beta = 1.0 / temperature;
  // calculate Boltzmann factor
  return std::pow(volume, current_reaction.nu_bar) * current_reaction.gamma *
         factorial_expr * exp(-beta * (E_pot_new - E_pot_old));
}

std::map<int, int>
ReactionAlgorithm::save_old_particle_numbers(int reaction_id) {
  std::map<int, int> old_particle_numbers;
  // reactants
  for (int type : reactions[reaction_id].reactant_types) {
    old_particle_numbers[type] = number_of_particles_with_type(type);
  }

  // products
  for (int type : reactions[reaction_id].product_types) {
    old_particle_numbers[type] = number_of_particles_with_type(type);
  }
  return old_particle_numbers;
}

void WangLandauReactionEnsemble::on_reaction_entry(int &old_state_index) {
  old_state_index = get_flattened_index_wang_landau_of_current_state();
  if (old_state_index >= 0) {
    if (histogram[old_state_index] >= 0)
      monte_carlo_trial_moves += 1;
  }
}

void WangLandauReactionEnsemble::on_reaction_rejection_directly_after_entry(
    int &old_state_index) {
  // increase the Wang-Landau potential and histogram at the current nbar
  // (this case covers the cases nbar=0 or nbar=1)
  update_wang_landau_potential_and_histogram(old_state_index);
}

void WangLandauReactionEnsemble::on_attempted_reaction(int &new_state_index) {
  new_state_index = get_flattened_index_wang_landau_of_current_state();
}

void WangLandauReactionEnsemble::on_end_reaction(int &accepted_state) {
  update_wang_landau_potential_and_histogram(accepted_state);
}

/**
 * Generic one way reaction
 * A+B+...+G +... --> K+...X + Z +...
 * You need to use 2A --> B instead of A+A --> B since in the last case you
 * assume distinctness of the particles, however both ways to describe the
 * reaction are equivalent in the thermodynamic limit (large particle numbers).
 * Furthermore, it is crucial for the function in which order you provide the
 * reactant and product types since particles will be replaced correspondingly!
 * If there are less reactants than products, new product particles are created
 * randomly in the box. Matching particles simply change the types. If there
 * are more reactants than products, old reactant particles are deleted.
 */
void ReactionAlgorithm::generic_oneway_reaction(int reaction_id) {

  SingleReaction &current_reaction = reactions[reaction_id];
  current_reaction.tried_moves += 1;
  particle_inside_exclusion_radius_touched = false;
  int old_state_index = -1; // for Wang-Landau algorithm
  on_reaction_entry(old_state_index);
  if (!all_reactant_particles_exist(reaction_id)) {
    // makes sure, no incomplete reaction is performed -> only need to consider
    // rollback of complete reactions
    on_reaction_rejection_directly_after_entry(old_state_index);
    return;
  }

  // calculate potential energy
  const double E_pot_old =
      calculate_current_potential_energy_of_system(); // only consider potential
                                                      // energy since we assume
                                                      // that the kinetic part
                                                      // drops out in the
                                                      // process of calculating
                                                      // ensemble averages
                                                      // (kinetic part may be
                                                      // separated and crossed
                                                      // out)

  // find reacting molecules in reactants and save their properties for later
  // recreation if step is not accepted
  // do reaction
  // save old particle_numbers
  std::map<int, int> old_particle_numbers =
      save_old_particle_numbers(reaction_id);

  std::vector<int> p_ids_created_particles;
  std::vector<StoredParticleProperty> hidden_particles_properties;
  std::vector<StoredParticleProperty> changed_particles_properties;
  const int number_of_saved_properties =
      3; // save p_id, charge and type of the reactant particle, only thing we
         // need to hide the particle and recover it
  make_reaction_attempt(current_reaction, changed_particles_properties,
                        p_ids_created_particles, hidden_particles_properties);

  double E_pot_new;
  if (particle_inside_exclusion_radius_touched)
    E_pot_new = std::numeric_limits<double>::max();
  else
    E_pot_new = calculate_current_potential_energy_of_system();

  int new_state_index = -1; // save new_state_index for Wang-Landau algorithm
  int accepted_state = -1;  // for Wang-Landau algorithm
  on_attempted_reaction(new_state_index);

  constexpr bool only_make_configuration_changing_move = false;
  auto const bf = calculate_acceptance_probability(
      current_reaction, E_pot_old, E_pot_new, old_particle_numbers,
      old_state_index, new_state_index, only_make_configuration_changing_move);

  std::vector<double> exponential = {
      exp(-1.0 / temperature * (E_pot_new - E_pot_old))};
  current_reaction.accumulator_exponentials(exponential);

  if (m_uniform_real_distribution(m_generator) < bf) {
    // accept
    accepted_state = new_state_index;

    // delete hidden reactant_particles (remark: don't delete changed particles)
    // extract ids of to be deleted particles
    auto len_hidden_particles_properties =
        static_cast<int>(hidden_particles_properties.size());
    std::vector<int> to_be_deleted_hidden_ids(len_hidden_particles_properties);
    std::vector<int> to_be_deleted_hidden_types(
        len_hidden_particles_properties);
    for (int i = 0; i < len_hidden_particles_properties; i++) {
      auto p_id = hidden_particles_properties[i].p_id;
      to_be_deleted_hidden_ids[i] = p_id;
      to_be_deleted_hidden_types[i] = hidden_particles_properties[i].type;
      // change back type otherwise the bookkeeping algorithm is not working
      set_particle_type(p_id, hidden_particles_properties[i].type);
    }

    for (int i = 0; i < len_hidden_particles_properties; i++) {
      delete_particle(to_be_deleted_hidden_ids[i]); // delete particle
    }
    current_reaction.accepted_moves += 1;
  } else {
    // reject
    accepted_state = old_state_index;
    // reverse reaction
    // 1) delete created product particles
    for (int p_ids_created_particle : p_ids_created_particles) {
      delete_particle(p_ids_created_particle);
    }
    // 2) restore previously hidden reactant particles
    restore_properties(hidden_particles_properties, number_of_saved_properties);
    // 3) restore previously changed reactant particles
    restore_properties(changed_particles_properties,
                       number_of_saved_properties);
  }
  on_end_reaction(accepted_state);
}

/**
 * Replaces a particle with the given particle id to be of a certain type. This
 * especially means that the particle type and the particle charge are changed.
 */
void ReactionAlgorithm::replace_particle(int p_id, int desired_type) {
  set_particle_type(p_id, desired_type);
#ifdef ELECTROSTATICS
  set_particle_q(p_id, charges_of_types[desired_type]);
#endif
}

/**
 * Hides a particle from short ranged interactions and from the electrostatic
 * interaction. Additional hiding from interactions would need to be implemented
 * here.
 *
 * Removing the particle charge and changing its type to a non existing one
 * deactivates all interactions with other particles, as if the particle was
 * inexistent (currently only type-based interactions are switched off, as well
 * as the electrostatic interaction).
 * This function does not break bonds for simple reactions, as long as there
 * are no reactions like 2A -->B where one of the reacting A particles occurs
 * in the polymer (think of bond breakages if the monomer in the polymer gets
 * deleted in the reaction). This constraint is not of fundamental reason, but
 * there would be a need for a rule for such "collision" reactions (a reaction
 * like the one above).
 */
void ReactionAlgorithm::hide_particle(int p_id, int previous_type) {

  auto const part = get_particle_data(p_id);
  auto const d_min = distto(partCfg(), part.r.p, p_id);
  if (d_min < exclusion_radius)
    particle_inside_exclusion_radius_touched = true;

#ifdef ELECTROSTATICS
  // set charge
  set_particle_q(p_id, 0.0);
#endif
  // set type
  set_particle_type(p_id, non_interacting_type);
}

/**
 * Deletes the particle with the given p_id and stores the id if the deletion
 * created a hole in the particle id range. This method is intended to only
 * delete unbonded particles since bonds are coupled to ids. This is used to
 * avoid the id range becoming excessively huge.
 */
int ReactionAlgorithm::delete_particle(int p_id) {
  auto const old_max_seen_id = get_maximal_particle_id();
  if (p_id == old_max_seen_id) {
    // last particle, just delete
    remove_particle(p_id);
    // remove all saved empty p_ids which are greater than the max_seen_particle
    // this is needed in order to avoid the creation of holes
    for (auto p_id_iter = m_empty_p_ids_smaller_than_max_seen_particle.begin();
         p_id_iter != m_empty_p_ids_smaller_than_max_seen_particle.end();) {
      if ((*p_id_iter) >= old_max_seen_id)
        p_id_iter = m_empty_p_ids_smaller_than_max_seen_particle.erase(
            p_id_iter); // update iterator after container was modified
      else
        ++p_id_iter;
    }
  } else if (p_id <= old_max_seen_id) {
    remove_particle(p_id);
    m_empty_p_ids_smaller_than_max_seen_particle.push_back(p_id);
  } else {
    throw std::runtime_error(
        "Particle id is greater than the max seen particle id");
  }
  return 0;
}

/**
 * Writes a random position inside the central box into the provided array.
 */
Utils::Vector3d ReactionAlgorithm::get_random_position_in_box() {
  Utils::Vector3d out_pos{};

  if (box_is_cylindric_around_z_axis) {
    // see http://mathworld.wolfram.com/DiskPointPicking.html
    double random_radius =
        cyl_radius *
        std::sqrt(m_uniform_real_distribution(
            m_generator)); // for uniform disk point picking in cylinder
    double phi = 2.0 * Utils::pi() * m_uniform_real_distribution(m_generator);
    out_pos[0] = cyl_x + random_radius * cos(phi);
    out_pos[1] = cyl_y + random_radius * sin(phi);
    out_pos[2] = box_geo.length()[2] * m_uniform_real_distribution(m_generator);
  } else if (box_has_wall_constraints) {
    out_pos[0] = box_geo.length()[0] * m_uniform_real_distribution(m_generator);
    out_pos[1] = box_geo.length()[1] * m_uniform_real_distribution(m_generator);
    out_pos[2] = slab_start_z + (slab_end_z - slab_start_z) *
                                    m_uniform_real_distribution(m_generator);
  } else {
    // cubic case
    out_pos[0] = box_geo.length()[0] * m_uniform_real_distribution(m_generator);
    out_pos[1] = box_geo.length()[1] * m_uniform_real_distribution(m_generator);
    out_pos[2] = box_geo.length()[2] * m_uniform_real_distribution(m_generator);
  }
  return out_pos;
}

/**
 * Creates a particle at the end of the observed particle id range.
 */
int ReactionAlgorithm::create_particle(int desired_type) {
  int p_id;
  if (!m_empty_p_ids_smaller_than_max_seen_particle.empty()) {
    auto p_id_iter = std::min_element(
        std::begin(m_empty_p_ids_smaller_than_max_seen_particle),
        std::end(m_empty_p_ids_smaller_than_max_seen_particle));
    p_id = *p_id_iter;
    m_empty_p_ids_smaller_than_max_seen_particle.erase(p_id_iter);
  } else {
    p_id = get_maximal_particle_id() + 1;
  }
  Utils::Vector3d pos_vec;

  // create random velocity vector according to Maxwell-Boltzmann distribution
  double vel[3];
  // we use mass=1 for all particles, think about adapting this
  vel[0] = std::sqrt(temperature) * m_normal_distribution(m_generator);
  vel[1] = std::sqrt(temperature) * m_normal_distribution(m_generator);
  vel[2] = std::sqrt(temperature) * m_normal_distribution(m_generator);
#ifdef ELECTROSTATICS
  double charge = charges_of_types[desired_type];
#endif

  pos_vec = get_random_position_in_box();
  place_particle(p_id, pos_vec.data());
  // set type
  set_particle_type(p_id, desired_type);
#ifdef ELECTROSTATICS
  // set charge
  set_particle_q(p_id, charge);
#endif
  // set velocities
  set_particle_v(p_id, vel);
  double d_min = distto(partCfg(), pos_vec, p_id);
  if (d_min < exclusion_radius) {
    // setting of a minimal distance is allowed to avoid overlapping
    // configurations if there is a repulsive potential. States with
    // very high energies have a probability of almost zero and
    // therefore do not contribute to ensemble averages.
    particle_inside_exclusion_radius_touched = true;
  }
  return p_id;
}

void WangLandauReactionEnsemble::on_mc_rejection_directly_after_entry(
    int &old_state_index) {
  if (do_energy_reweighting)
    update_wang_landau_potential_and_histogram(old_state_index);
}

void WangLandauReactionEnsemble::on_mc_accept(int &new_state_index) {
  if (do_energy_reweighting) {
    // modify wang_landau histogram and potential
    update_wang_landau_potential_and_histogram(new_state_index);
  }
}

void WangLandauReactionEnsemble::on_mc_reject(int &old_state_index) {
  if (do_energy_reweighting)
    update_wang_landau_potential_and_histogram(old_state_index);
}

int WangLandauReactionEnsemble::on_mc_use_WL_get_new_state() {
  return get_flattened_index_wang_landau_of_current_state();
}

/**
 * Performs a global MC move for a particle of the provided type.
 */
bool ReactionAlgorithm::do_global_mc_move_for_particles_of_type(
    int type, int particle_number_of_type_to_be_changed, bool use_wang_landau) {
  m_tried_configurational_MC_moves += 1;
  particle_inside_exclusion_radius_touched = false;

  int old_state_index = -1;
  if (use_wang_landau) {
    on_reaction_entry(old_state_index);
  }

  int particle_number_of_type = number_of_particles_with_type(type);
  if (particle_number_of_type == 0 or
      particle_number_of_type_to_be_changed == 0) {
    // reject
    if (use_wang_landau) {
      on_mc_rejection_directly_after_entry(old_state_index);
    }
    return false;
  }

  const double E_pot_old = calculate_current_potential_energy_of_system();

  std::vector<double> particle_positions(3 *
                                         particle_number_of_type_to_be_changed);
  std::vector<int> p_id_s_changed_particles;

  // heuristic to find a particle that hasn't been touched already
  int random_index_in_type_map = i_random(number_of_particles_with_type(type));
  int p_id = get_random_p_id(type, random_index_in_type_map);
  for (int i = 0; i < particle_number_of_type_to_be_changed; i++) {
    while (Utils::contains(p_id_s_changed_particles, p_id)) {
      random_index_in_type_map = i_random(number_of_particles_with_type(type));
      p_id = get_random_p_id(
          type, random_index_in_type_map); // check whether you already touched
                                           // this p_id, then reassign
    }

    auto part = get_particle_data(p_id);

    particle_positions[3 * i] = part.r.p[0];
    particle_positions[3 * i + 1] = part.r.p[1];
    particle_positions[3 * i + 2] = part.r.p[2];
    p_id_s_changed_particles.push_back(p_id);
  }

  // propose new positions
  for (int i = 0; i < particle_number_of_type_to_be_changed; i++) {
    p_id = p_id_s_changed_particles[i];
    // change particle position
    auto const new_pos = get_random_position_in_box();
    auto const &p = get_particle_data(p_id);
    auto const prefactor = std::sqrt(temperature / p.p.mass);
    double vel[3];
    vel[0] = prefactor * m_normal_distribution(m_generator);
    vel[1] = prefactor * m_normal_distribution(m_generator);
    vel[2] = prefactor * m_normal_distribution(m_generator);
    set_particle_v(p_id, vel);
    place_particle(p_id, new_pos.data());
    auto const d_min = distto(partCfg(), new_pos, p_id);
    if (d_min < exclusion_radius)
      particle_inside_exclusion_radius_touched = true;
  }

  double E_pot_new;
  if (particle_inside_exclusion_radius_touched)
    E_pot_new = std::numeric_limits<double>::max();
  else
    E_pot_new = calculate_current_potential_energy_of_system();

  double beta = 1.0 / temperature;

  int new_state_index = -1;
  double bf = 1.0;
  std::map<int, int> dummy_old_particle_numbers;
  SingleReaction temp_unimportant_arbitrary_reaction;

  if (use_wang_landau) {
    new_state_index = on_mc_use_WL_get_new_state();
    bf = calculate_acceptance_probability(
        temp_unimportant_arbitrary_reaction, E_pot_old, E_pot_new,
        dummy_old_particle_numbers, old_state_index, new_state_index, true);
  } else {
    // Metropolis algorithm since proposal density is symmetric
    bf = std::min(1.0, bf * exp(-beta * (E_pot_new - E_pot_old)));
  }

  // // correct for enhanced proposal of small radii by using the
  // // Metropolis-Hastings algorithm for asymmetric proposal densities
  // double old_radius =
  //     std::sqrt(std::pow(particle_positions[0]-cyl_x,2) +
  //               std::pow(particle_positions[1]-cyl_y,2));
  // double new_radius =
  //     std::sqrt(std::pow(new_pos[0]-cyl_x,2)+std::pow(new_pos[1]-cyl_y,2));
  // bf = std::min(1.0,
  //     bf*exp(-beta*(E_pot_new-E_pot_old))*new_radius/old_radius);

  // Metropolis-Hastings algorithm for asymmetric proposal density
  if (m_uniform_real_distribution(m_generator) < bf) {
    // accept
    m_accepted_configurational_MC_moves += 1;
    if (use_wang_landau) {
      on_mc_accept(new_state_index);
    }
    return true;
  }
  // reject
  // modify wang_landau histogram and potential
  if (use_wang_landau) {
    on_mc_reject(old_state_index);
  }
  // create particles again at the positions they were
  for (int i = 0; i < particle_number_of_type_to_be_changed; i++)
    place_particle(p_id_s_changed_particles[i], &particle_positions[3 * i]);
  return false;
}

/**
 * Adds a new collective variable (CV) of the type degree of association to the
 * Wang-Landau sampling
 */
void WangLandauReactionEnsemble::add_new_CV_degree_of_association(
    int associated_type, double CV_minimum, double CV_maximum,
    const std::vector<int> &corresponding_acid_types) {
  std::shared_ptr<DegreeOfAssociationCollectiveVariable>
      new_collective_variable =
          std::make_shared<DegreeOfAssociationCollectiveVariable>();
  new_collective_variable->associated_type = associated_type;
  new_collective_variable->CV_minimum = CV_minimum;
  new_collective_variable->CV_maximum = CV_maximum;
  new_collective_variable->corresponding_acid_types = corresponding_acid_types;
  new_collective_variable->delta_CV =
      calculate_delta_degree_of_association(*new_collective_variable);
  collective_variables.push_back(new_collective_variable);
  initialize_wang_landau();
}

/**
 * Adds a new collective variable (CV) of the type potential energy to the
 * Wang-Landau sampling
 */
void WangLandauReactionEnsemble::add_new_CV_potential_energy(
    std::istream &infile, double delta_CV) {
  std::shared_ptr<EnergyCollectiveVariable> new_collective_variable =
      std::make_shared<EnergyCollectiveVariable>();
  new_collective_variable->delta_CV = delta_CV;
  new_collective_variable->load_CV_boundaries(*this, infile);
  collective_variables.emplace_back(new_collective_variable);
  initialize_wang_landau();
}

/**
 * Adds a new collective variable (CV) of the type potential energy to the
 * Wang-Landau sampling
 */
void WangLandauReactionEnsemble::add_new_CV_potential_energy(
    const std::string &filename, double delta_CV) {
  std::ifstream infile;
  infile.open(filename);
  if (!infile.is_open()) {
    throw std::runtime_error("Cannot read " + filename);
  }
  add_new_CV_potential_energy(infile, delta_CV);
  infile.close();
}

/**
 * Returns the flattened index of the multidimensional Wang-Landau histogram
 */
int WangLandauReactionEnsemble::get_flattened_index_wang_landau(
    std::vector<double> const &current_state,
    std::vector<double> const &collective_variables_minimum_values,
    std::vector<double> const &collective_variables_maximum_values,
    std::vector<double> const &delta_collective_variables_values,
    int nr_collective_variables) {

  // check for the current state to be an allowed state in the range
  // [collective_variables_minimum_values:collective_variables_maximum_values],
  // else return a negative index (to indicate error)
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    if (current_state[CV_i] >
            collective_variables_maximum_values[CV_i] +
                delta_collective_variables_values[CV_i] * 0.98 ||
        current_state[CV_i] <
            collective_variables_minimum_values[CV_i] -
                delta_collective_variables_values[CV_i] * 0.01) {
      return -10;
    }
  }

  std::vector<int> individual_indices(nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    auto const value =
        (current_state[CV_i] - collective_variables_minimum_values[CV_i]) /
        delta_collective_variables_values[CV_i];
    int rounded_value;
    if (CV_i == collective_variables.size() - 1 && do_energy_reweighting) {
      // for energy collective variable (simple truncating conversion desired)
      rounded_value = static_cast<int>(value);
    } else {
      // for degree of association collective variables (rounding conversion)
      rounded_value = static_cast<int>(std::floor(value));
    }
    if (rounded_value < 0 or
        rounded_value >= nr_subindices_of_collective_variable[CV_i]) {
      return -10;
    }
    individual_indices[CV_i] = rounded_value;
  }
  // get flattened index from individual_indices
  int index = 0;
  // this is already part of the algorithm to find the correct index
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    int factor = 1;
    for (int j = CV_i + 1; j < nr_collective_variables; j++) {
      factor *= nr_subindices_of_collective_variable[j];
    }
    index += factor * individual_indices[CV_i];
  }
  return index;
}

/**
 * Returns the flattened index of the multidimensional Wang-Landau histogram for
 * the current state of the simulation.
 */
int WangLandauReactionEnsemble::
    get_flattened_index_wang_landau_of_current_state() {
  auto const nr_collective_variables =
      static_cast<int>(collective_variables.size());
  // get current state
  std::vector<double> current_state(nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++)
    current_state[CV_i] = collective_variables[CV_i]->determine_current_state();
  // get collective_variables_minimum_values
  std::vector<double> collective_variables_minimum_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    collective_variables_minimum_values[CV_i] =
        collective_variables[CV_i]->CV_minimum;
  }
  // get collective_variables_maximum_values
  std::vector<double> collective_variables_maximum_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    collective_variables_maximum_values[CV_i] =
        collective_variables[CV_i]->CV_maximum;
  }
  // get delta_collective_variables_values
  std::vector<double> delta_collective_variables_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    delta_collective_variables_values[CV_i] =
        collective_variables[CV_i]->delta_CV;
  }
  int index = get_flattened_index_wang_landau(
      current_state, collective_variables_minimum_values,
      collective_variables_maximum_values, delta_collective_variables_values,
      nr_collective_variables);
  return index;
}

/**
 * Returns the minimum value of the collective variable on a delta_CV spaced
 * grid which starts at 0
 */
double WangLandauReactionEnsemble::get_minimum_CV_value_on_delta_CV_spaced_grid(
    double min_CV_value, double delta_CV) const {
  // assume grid has it s origin at 0
  double minimum_CV_value_on_delta_CV_spaced_grid =
      floor(min_CV_value / delta_CV) * delta_CV;
  return minimum_CV_value_on_delta_CV_spaced_grid;
}

/**
 * Calculates the smallest difference in the degree of association which can be
 * observed when changing the degree of association by one single reaction.
 */
double WangLandauReactionEnsemble::calculate_delta_degree_of_association(
    DegreeOfAssociationCollectiveVariable &current_collective_variable) {
  // calculate Delta in the degree of association so that EVERY reaction step is
  // driven.
  int total_number_of_corresponding_acid = 0;
  for (int corresponding_acid_type :
       current_collective_variable.corresponding_acid_types) {
    int num_of_current_type =
        number_of_particles_with_type(corresponding_acid_type);
    total_number_of_corresponding_acid += num_of_current_type;
  }
  double delta = 1.0 / total_number_of_corresponding_acid;
  // now modify the minimum value of the CV to lie on the grid
  current_collective_variable.CV_minimum =
      get_minimum_CV_value_on_delta_CV_spaced_grid(
          current_collective_variable.CV_minimum, delta);
  return delta;
}

void WangLandauReactionEnsemble::invalidate_bins() {
  // make values in histogram and Wang-Landau potential negative if they are not
  // allowed at the given degree of association, because the energy boundaries
  // prohibit them

  long empty_bins_in_memory = 0;
  for (std::size_t flattened_index = 0;
       flattened_index < wang_landau_potential.size(); flattened_index++) {
    std::vector<std::size_t> unraveled_index(
        nr_subindices_of_collective_variable.size());
    Utils::unravel_index(nr_subindices_of_collective_variable.begin(),
                         nr_subindices_of_collective_variable.end(),
                         unraveled_index.begin(), flattened_index);
    // use unraveled index
    int EnergyCollectiveVariable_index = 0;
    if (collective_variables.size() > 1) {
      // assume the energy collective variable to be the last added
      // collective variable
      EnergyCollectiveVariable_index =
          static_cast<int>(collective_variables.size()) - 1;
    }
    double current_energy =
        static_cast<double>(unraveled_index[EnergyCollectiveVariable_index]) *
            collective_variables[EnergyCollectiveVariable_index]->delta_CV +
        collective_variables[EnergyCollectiveVariable_index]->CV_minimum;
    int flat_index_without_energy_CV =
        get_flattened_index_wang_landau_without_energy_collective_variable(
            static_cast<int>(flattened_index), EnergyCollectiveVariable_index);
    std::shared_ptr<CollectiveVariable> energy_CV =
        collective_variables[EnergyCollectiveVariable_index];
    if (current_energy >
            max_boundaries_energies[flat_index_without_energy_CV] ||
        current_energy < min_boundaries_energies[flat_index_without_energy_CV] -
                             energy_CV->delta_CV) {
      histogram[flattened_index] = int_fill_value;
      wang_landau_potential[flattened_index] = double_fill_value;
      empty_bins_in_memory++;
    }
  }

  used_bins =
      static_cast<long>(wang_landau_potential.size()) - empty_bins_in_memory;
}

int WangLandauReactionEnsemble::initialize_wang_landau() {

  nr_subindices_of_collective_variable.resize(collective_variables.size(), 0);
  auto const &cv = collective_variables.back();
  // add 1 for collective variables which are of type degree of association
  nr_subindices_of_collective_variable.back() =
      static_cast<int>((cv->CV_maximum - cv->CV_minimum) / cv->delta_CV) + 1;

  // construct (possibly higher dimensional) histogram and potential over
  // gamma (the space which should be equally sampled when the Wang-Landau
  // algorithm has converged)
  // add 1 for degrees of association-related part of histogram (think of only
  // one acid particle)
  auto const needed_bins = std::accumulate(
      collective_variables.begin(), collective_variables.end(), 1l,
      [](long acc, std::shared_ptr<CollectiveVariable> const &cv) {
        return acc * (static_cast<long>((cv->CV_maximum - cv->CV_minimum) /
                                        cv->delta_CV) +
                      1l);
      });
  assert(needed_bins >= 0);
  histogram.resize(needed_bins, 0);
  wang_landau_potential.resize(needed_bins, 0);

  used_bins = needed_bins; // initialize for 1/t wang_landau algorithm

  if (do_energy_reweighting) {
    invalidate_bins();
  }
  return ES_OK;
}

/** Calculate the expression in the acceptance probability of the Wang-Landau
 *  reaction ensemble.
 *  Modify Boltzmann factor according to Wang-Landau algorithm in @cite yan02b.
 */
double WangLandauReactionEnsemble::calculate_acceptance_probability(
    SingleReaction const &current_reaction, double E_pot_old, double E_pot_new,
    std::map<int, int> const &old_particle_numbers, int old_state_index,
    int new_state_index, bool only_make_configuration_changing_move) const {

  double beta = 1.0 / temperature;
  double bf = 1.0;

  if (!(do_not_sample_reaction_partition_function ||
        only_make_configuration_changing_move)) {
    auto const factorial_expr =
        calculate_factorial_expression(current_reaction, old_particle_numbers);
    bf = std::pow(volume, current_reaction.nu_bar) * current_reaction.gamma *
         factorial_expr;
  }

  if (!do_energy_reweighting) {
    bf *= exp(-beta * (E_pot_new - E_pot_old));
  }

  // Check whether the proposed state lies in the reaction coordinate space
  // gamma and add the Wang-Landau modification factor, this is a bit nasty
  // due to the energy collective variable case (memory layout of storage
  // array of the histogram and the wang_landau_potential values is "cuboid").
  if (old_state_index >= 0 && new_state_index >= 0) {
    if (histogram[new_state_index] >= 0 && histogram[old_state_index] >= 0) {
      // Modify Boltzmann factor (bf) according to Wang-Landau algorithm
      // @cite yan02b.
      // This makes the new state being accepted with the conditional
      // probability bf (bf is a transition probability = conditional
      // probability from the old state to move to the new state).
      bf = std::min(1.0, bf * exp(wang_landau_potential[old_state_index] -
                                  wang_landau_potential[new_state_index]));
    } else {
      if (histogram[old_state_index] < 0)
        bf = 10; // accept the reaction if we found a state in gamma
                 // (histogram[new_state_index] >= 0) or to sample new configs
                 // which might lie in gamma(histogram[new_state_index] < 0)
      else if (histogram[new_state_index] < 0)
        bf = -10; // reject the reaction, since the new state
                  // is not in gamma while the old sate was in gamma
    }
  } else if (old_state_index < 0 && new_state_index >= 0) {
    bf = 10; // accept the reaction if we found a new state in gamma
             // (new_state_index >= 0) or to sample new configs which
             // might lie in gamma (new_state_index < 0)
  } else if (old_state_index >= 0 && new_state_index < 0) {
    bf = -10; // reject the reaction, since the new state is
              // not in gamma while the old sate was in gamma
  }
  return bf;
}

/** Perform a randomly selected reaction using the Wang-Landau algorithm.
 *
 *  Make sure to perform additional configuration changing steps, after the
 *  reaction step! Like in @cite yan02b. This can be done with MD in the case
 *  of the no-energy-reweighting case, or with the function
 *  @ref ReactionAlgorithm::do_global_mc_move_for_particles_of_type.
 *
 *  Perform additional Monte Carlo moves to sample configurational
 *  partition function according to @cite yan02b. Do as many steps
 *  as needed to get to a new conformation.
 */
int WangLandauReactionEnsemble::do_reaction(int reaction_steps) {
  m_WL_tries += reaction_steps;
  for (int step = 0; step < reaction_steps; step++) {
    int reaction_id = i_random(static_cast<int>(reactions.size()));
    generic_oneway_reaction(reaction_id);
    if (can_refine_wang_landau_one_over_t() && m_WL_tries % 10000 == 0) {
      // check for convergence
      if (achieved_desired_number_of_refinements_one_over_t()) {
        write_wang_landau_results_to_file(output_filename);
        return -10; // return negative value to indicate that the Wang-Landau
                    // algorithm has converged
      }
      refine_wang_landau_parameter_one_over_t();
    }
  }
  // shift Wang-Landau potential minimum to zero
  if (m_WL_tries % (std::max(90000, 9 * reaction_steps)) == 0) {
    // for numerical stability here we also subtract the minimum positive value
    // of the wang_landau_potential from the wang_landau potential, allowed
    // since only the difference in the wang_landau potential is of interest.
    double minimum_wang_landau_potential =
        find_minimum_non_negative_value(wang_landau_potential);
    for (double &i : wang_landau_potential) {
      if (i >= 0) // check for whether we are in the
                  // valid range of the collective variable
        i -= minimum_wang_landau_potential;
    }
    // write out preliminary Wang-Landau potential results
    write_wang_landau_results_to_file(output_filename);
  }
  return 0;
}

/** Increase the Wang-Landau potential and histogram at the current nbar */
void WangLandauReactionEnsemble::update_wang_landau_potential_and_histogram(
    int index_of_state_after_acceptance_or_rejection) {
  if (index_of_state_after_acceptance_or_rejection >= 0) {
    if (histogram[index_of_state_after_acceptance_or_rejection] >= 0) {
      histogram[index_of_state_after_acceptance_or_rejection] += 1;
      wang_landau_potential[index_of_state_after_acceptance_or_rejection] +=
          wang_landau_parameter;
    }
  }
}

/**
 *Determines whether we can reduce the Wang-Landau parameter
 */
bool WangLandauReactionEnsemble::can_refine_wang_landau_one_over_t() const {
  double minimum_required_value =
      0.80 * average_list_of_allowed_entries(
                 histogram); // This is an additional constraint to sample
                             // configuration space better. Use flatness
                             // criterion according to 1/t algorithm as long as
                             // you are not in 1/t regime.
  if (do_energy_reweighting)
    minimum_required_value = 20; // get faster in energy reweighting case

  return *(std::min_element(histogram.begin(), histogram.end())) >
             minimum_required_value ||
         m_system_is_in_1_over_t_regime;
}

/**
 *Reset the Wang-Landau histogram.
 */
void WangLandauReactionEnsemble::reset_histogram() {
  printf("Histogram is flat. Refining. Previous Wang-Landau modification "
         "parameter was %f.\n",
         wang_landau_parameter);
  fflush(stdout);

  for (int i = 0; i < wang_landau_potential.size(); i++) {
    if (histogram[i] >= 0) { // checks for validity of index i (think of energy
                             // collective variables, in a cubic memory layout
                             // there will be indices which are not allowed by
                             // the energy boundaries. These values will be
                             // initialized with a negative fill value)
      histogram[i] = 0;
    }
  }
}

/**
 *Refine the Wang-Landau parameter using the 1/t rule.
 */
void WangLandauReactionEnsemble::refine_wang_landau_parameter_one_over_t() {
  double monte_carlo_time = static_cast<double>(monte_carlo_trial_moves) /
                            static_cast<double>(used_bins);
  if (wang_landau_parameter / 2.0 <= 1.0 / monte_carlo_time ||
      m_system_is_in_1_over_t_regime) {
    wang_landau_parameter = 1.0 / monte_carlo_time;
    if (!m_system_is_in_1_over_t_regime) {
      m_system_is_in_1_over_t_regime = true;
      printf("Refining: Wang-Landau parameter is now 1/t.\n");
    }
  } else {
    reset_histogram();
    wang_landau_parameter = wang_landau_parameter / 2.0;
  }
}

/**
 *Determine whether the desired number of refinements was achieved.
 */
bool WangLandauReactionEnsemble::
    achieved_desired_number_of_refinements_one_over_t() const {
  if (wang_landau_parameter < final_wang_landau_parameter) {
    printf("Achieved desired number of refinements\n");
    return true;
  }
  return false;
}

void WangLandauReactionEnsemble::format_wang_landau_results(std::ostream &out) {
  for (std::size_t flattened_index = 0;
       flattened_index < wang_landau_potential.size(); flattened_index++) {
    // unravel index
    if (std::abs(wang_landau_potential[flattened_index] - double_fill_value) >
        1) {
      // only output data if they are not equal to double_fill_value.
      // This ensures that for the energy observable not allowed energies
      // (energies in the interval [global_E_min, global_E_max]) in the
      // multidimensional Wang-Landau potential are printed out, since
      // the range [E_min(nbar), E_max(nbar)] for each nbar may be a
      // different one
      std::vector<std::size_t> unraveled_index(
          nr_subindices_of_collective_variable.size());
      Utils::unravel_index(nr_subindices_of_collective_variable.begin(),
                           nr_subindices_of_collective_variable.end(),
                           unraveled_index.begin(), flattened_index);
      // use unraveled index
      for (std::size_t i = 0; i < collective_variables.size(); i++) {
        auto const value = static_cast<double>(unraveled_index[i]) *
                               collective_variables[i]->delta_CV +
                           collective_variables[i]->CV_minimum;
        out << value << " ";
      }
      out << wang_landau_potential[flattened_index] << " \n";
    }
  }
}

/**
 *Writes the Wang-Landau potential to file.
 */
void WangLandauReactionEnsemble::write_wang_landau_results_to_file(
    const std::string &filename) {
  std::ofstream outfile;

  outfile.open(filename);
  if (!outfile.is_open()) {
    throw std::runtime_error("Cannot write to " + filename);
  }

  format_wang_landau_results(outfile);
  outfile.close();
}

/**
 *Update the minimum and maximum observed energies using the current state.
 *Needed for preliminary energy reweighting runs.
 */
int WangLandauReactionEnsemble::
    update_maximum_and_minimum_energies_at_current_state() {
  if (minimum_energies_at_flat_index.empty() ||
      maximum_energies_at_flat_index.empty()) {
    minimum_energies_at_flat_index.resize(wang_landau_potential.size(),
                                          double_fill_value);
    maximum_energies_at_flat_index.resize(wang_landau_potential.size(),
                                          double_fill_value);
  }

  const double E_pot_current = calculate_current_potential_energy_of_system();
  int index = get_flattened_index_wang_landau_of_current_state();

  // update stored energy values
  if (((E_pot_current < minimum_energies_at_flat_index[index]) ||
       std::abs(minimum_energies_at_flat_index[index] - double_fill_value) <
           std::numeric_limits<double>::epsilon())) {
    minimum_energies_at_flat_index[index] = E_pot_current;
  }
  if (((E_pot_current > maximum_energies_at_flat_index[index]) ||
       std::abs(maximum_energies_at_flat_index[index] - double_fill_value) <
           std::numeric_limits<double>::epsilon())) {
    maximum_energies_at_flat_index[index] = E_pot_current;
  }

  return 0;
}

/**
 *Write out an energy boundary file using the energy boundaries observed in a
 *preliminary energy reweighting run.
 */
void WangLandauReactionEnsemble::write_out_preliminary_energy_run_results(
    const std::string &filename) {

  std::ofstream outfile;
  outfile.open(filename);
  if (!outfile.is_open()) {
    throw std::runtime_error("Cannot write to " + filename);
  }
  outfile << "#nbar E_min E_max\n";
  for (std::size_t flattened_index = 0;
       flattened_index < wang_landau_potential.size(); flattened_index++) {
    // unravel index
    std::vector<std::size_t> unraveled_index(
        nr_subindices_of_collective_variable.size());
    Utils::unravel_index(nr_subindices_of_collective_variable.begin(),
                         nr_subindices_of_collective_variable.end(),
                         unraveled_index.begin(), flattened_index);
    // use unraveled index
    for (std::size_t i = 0; i < collective_variables.size(); i++) {
      auto const value = static_cast<double>(unraveled_index[i]) *
                             collective_variables[i]->delta_CV +
                         collective_variables[i]->CV_minimum;
      outfile << value << " ";
    }
    outfile << minimum_energies_at_flat_index[flattened_index] << " "
            << maximum_energies_at_flat_index[flattened_index] << " \n";
  }
  outfile.close();
}

/**
 *Returns the flattened index of a given flattened index without the energy
 *collective variable.
 */
int WangLandauReactionEnsemble::
    get_flattened_index_wang_landau_without_energy_collective_variable(
        int flattened_index_with_EnergyCollectiveVariable,
        int CV_index_energy_observable) {
  std::vector<std::size_t> unraveled_index(
      nr_subindices_of_collective_variable.size());
  Utils::unravel_index(
      nr_subindices_of_collective_variable.begin(),
      nr_subindices_of_collective_variable.end(), unraveled_index.begin(),
      static_cast<std::size_t>(flattened_index_with_EnergyCollectiveVariable));
  // use unraveled index, but skip last collective variable
  // (the energy collective variable)
  auto const nr_collective_variables =
      static_cast<int>(collective_variables.size()) - 1;
  std::vector<double> current_state(nr_collective_variables);
  for (int i = 0; i < nr_collective_variables; i++) {
    current_state[i] = static_cast<double>(unraveled_index[i]) *
                           collective_variables[i]->delta_CV +
                       collective_variables[i]->CV_minimum;
  }

  // get collective_variables_minimum_values
  std::vector<double> collective_variables_minimum_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    collective_variables_minimum_values[CV_i] =
        collective_variables[CV_i]->CV_minimum;
  }
  // get collective_variables_maximum_values
  std::vector<double> collective_variables_maximum_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    collective_variables_maximum_values[CV_i] =
        collective_variables[CV_i]->CV_maximum;
  }
  // get delta_collective_variables_values
  std::vector<double> delta_collective_variables_values(
      nr_collective_variables);
  for (int CV_i = 0; CV_i < nr_collective_variables; CV_i++) {
    delta_collective_variables_values[CV_i] =
        collective_variables[CV_i]->delta_CV;
  }
  int index = get_flattened_index_wang_landau(
      current_state, collective_variables_minimum_values,
      collective_variables_maximum_values, delta_collective_variables_values,
      nr_collective_variables);
  return index;
}

/**
 * Writes the Wang-Landau parameter, the histogram and the potential to a file.
 * You can restart a Wang-Landau simulation using this information.
 * Additionally you should store the positions of the particles.
 * Not storing them introduces small, small statistical errors.
 */
void WangLandauReactionEnsemble::write_wang_landau_checkpoint(
    const std::string &identifier) {
  std::string filename;
  std::ofstream outfile;

  // write current Wang-Landau parameters (wang_landau_parameter,
  // monte_carlo_trial_moves, flat_index_of_current_state)
  filename = std::string("checkpoint_wang_landau_parameters_") + identifier;
  outfile.open(filename);
  if (!outfile.is_open()) {
    throw std::runtime_error("Cannot write to " + filename);
  }
  outfile << wang_landau_parameter << " " << monte_carlo_trial_moves << " "
          << get_flattened_index_wang_landau_of_current_state() << "\n";
  outfile.close();

  // write histogram
  filename = std::string("checkpoint_wang_landau_histogram_") + identifier;
  outfile.open(filename);
  if (!outfile.is_open()) {
    throw std::runtime_error("Cannot write to " + filename);
  }
  for (int i = 0; i < wang_landau_potential.size(); i++) {
    outfile << histogram[i] << "\n";
  }
  outfile.close();
  // write Wang-Landau potential
  filename = std::string("checkpoint_wang_landau_potential_") + identifier;
  outfile.open(filename);
  if (!outfile.is_open()) {
    throw std::runtime_error("Cannot write to " + filename);
  }
  for (double i : wang_landau_potential) {
    outfile << i << "\n";
  }
  outfile.close();
}

/**
 *Loads the Wang-Landau checkpoint
 */
void WangLandauReactionEnsemble::load_wang_landau_checkpoint(
    const std::string &identifier) {
  std::string filename;
  std::ifstream infile;

  // restore Wang-Landau parameters
  filename = std::string("checkpoint_wang_landau_parameters_") + identifier;
  infile.open(filename);
  if (infile.is_open()) {
    double wang_landau_parameter_entry;
    int wang_landau_monte_carlo_trial_moves_entry;
    int flat_index_of_state_at_checkpointing;
    int line = 0;
    while (infile >> wang_landau_parameter_entry >>
           wang_landau_monte_carlo_trial_moves_entry >>
           flat_index_of_state_at_checkpointing) {
      wang_landau_parameter = wang_landau_parameter_entry;
      monte_carlo_trial_moves = wang_landau_monte_carlo_trial_moves_entry;
      line += 1;
    }
    infile.close();
  } else {
    throw std::runtime_error("Cannot read " + filename);
  }

  // restore histogram
  filename = std::string("checkpoint_wang_landau_histogram_") + identifier;
  infile.open(filename);
  if (infile.is_open()) {
    int hist_entry;
    int line = 0;
    while (infile >> hist_entry) {
      histogram[line] = hist_entry;
      line += 1;
    }
    infile.close();
  } else {
    throw std::runtime_error("Cannot read " + filename);
  }

  // restore Wang-Landau potential
  filename = std::string("checkpoint_wang_landau_potential_") + identifier;
  infile.open(filename);
  if (infile.is_open()) {
    double wang_landau_potential_entry;
    int line = 0;
    while (infile >> wang_landau_potential_entry) {
      wang_landau_potential[line] = wang_landau_potential_entry;
      line += 1;
    }
    infile.close();
  } else {
    throw std::runtime_error("Cannot read " + filename);
  }

  // possible task: restore state in which the system was when the checkpoint
  // was written. However as long as checkpointing and restoring the system form
  // the checkpoint is rare this should not matter statistically.
}

int ConstantpHEnsemble::get_random_valid_p_id() {
  int random_p_id = i_random(get_maximal_particle_id() + 1);
  // draw random p_ids till we draw a pid which exists
  while (not particle_exists(random_p_id))
    random_p_id = i_random(get_maximal_particle_id() + 1);
  return random_p_id;
}

/**
 *Performs a reaction in the constant pH ensemble
 */
int ConstantpHEnsemble::do_reaction(int reaction_steps) {


  for (int reaction_id = 0; reaction_id < reactions.size(); reaction_id++) {
    if (all_reactant_particles_exist(reaction_id) == false){return 1;
    }
  }
      
  for (int i = 0; i < reaction_steps; ++i) {
    // get a list of reactions where a randomly selected particle type occurs in
    // the reactant list. the selection probability of the particle types has to
    // be proportional to the number of occurrences of the number of particles
    // with
    // this type

    // for optimizations this list could be determined during the initialization
    std::vector<int> list_of_reaction_ids_with_given_reactant_type;
    while (list_of_reaction_ids_with_given_reactant_type
               .empty()) { // avoid selecting a (e.g. salt) particle which
                           // does not take part in a reaction
      int random_p_id = get_random_valid_p_id(); // only used to determine which
                                                 // reaction is attempted.
      auto part = get_particle_data(random_p_id);

      int type_of_random_p_id = part.p.type;

      // construct list of reactions with the above reactant type
      for (int reaction_i = 0; reaction_i < reactions.size(); reaction_i++) {

        SingleReaction &current_reaction = reactions[reaction_i];
        for (int reactant_i = 0; reactant_i < 1;
             reactant_i++) { // reactant_i < 1 since it is assumed in this place
                             // that the types A and HA occur in the first place
                             // only. These are the types that should be
                             // switched, H+ should not be switched
          if (current_reaction.reactant_types[reactant_i] ==
              type_of_random_p_id) {
            list_of_reaction_ids_with_given_reactant_type.push_back(reaction_i);
            break;
          }
        }
      }
    }

    // randomly select a reaction to be performed
    int reaction_id =
        list_of_reaction_ids_with_given_reactant_type[i_random(static_cast<int>(
            list_of_reaction_ids_with_given_reactant_type.size()))];
    generic_oneway_reaction(reaction_id);
  }
  return 0;
}

/**
 * Calculates the expression in the acceptance probability of the constant pH
 * method.
 */
double ConstantpHEnsemble::calculate_acceptance_probability(
    SingleReaction const &current_reaction, double E_pot_old, double E_pot_new,
    std::map<int, int> const &dummy_old_particle_numbers,
    int dummy_old_state_index, int dummy_new_state_index,
    bool dummy_only_make_configuration_changing_move) const {
  auto const beta = 1.0 / temperature;
  auto const pKa = -current_reaction.nu_bar * log10(current_reaction.gamma);
  auto const ln_bf = (E_pot_new - E_pot_old) - current_reaction.nu_bar / beta *
                                                   log(10) *
                                                   (m_constant_pH - pKa);
  auto const bf = exp(-beta * ln_bf);
  return bf;
}

std::pair<double, double>
WidomInsertion::measure_excess_chemical_potential(int reaction_id) {
  if (!all_reactant_particles_exist(reaction_id))
    throw std::runtime_error("Trying to remove some non-existing particles "
                             "from the system via the inverse Widom scheme.");

  SingleReaction &current_reaction = reactions[reaction_id];
  const double E_pot_old = calculate_current_potential_energy_of_system();

  // make reaction attempt
  std::vector<int> p_ids_created_particles;
  std::vector<StoredParticleProperty> hidden_particles_properties;
  std::vector<StoredParticleProperty> changed_particles_properties;
  const int number_of_saved_properties =
      3; // save p_id, charge and type of the reactant particle, only thing we
         // need to hide the particle and recover it
  make_reaction_attempt(current_reaction, changed_particles_properties,
                        p_ids_created_particles, hidden_particles_properties);
  const double E_pot_new = calculate_current_potential_energy_of_system();
  // reverse reaction attempt
  // reverse reaction
  // 1) delete created product particles
  for (int p_ids_created_particle : p_ids_created_particles) {
    delete_particle(p_ids_created_particle);
  }
  // 2) restore previously hidden reactant particles
  restore_properties(hidden_particles_properties, number_of_saved_properties);
  // 3) restore previously changed reactant particles
  restore_properties(changed_particles_properties, number_of_saved_properties);
  std::vector<double> exponential = {
      exp(-1.0 / temperature * (E_pot_new - E_pot_old))};
  current_reaction.accumulator_exponentials(exponential);

  // calculate mean excess chemical potential and standard error of the mean
  std::pair<double, double> result = std::make_pair(
      -temperature * log(current_reaction.accumulator_exponentials.mean()[0]),
      std::abs(-temperature /
               current_reaction.accumulator_exponentials.mean()[0] *
               current_reaction.accumulator_exponentials.std_error()[0]));
  return result;
}

/**
 * Calculates the whole product of factorial expressions which occur in the
 * reaction ensemble acceptance probability.
 *
 * See @cite smith94c.
 */
double
calculate_factorial_expression(SingleReaction const &current_reaction,
                               std::map<int, int> const &old_particle_numbers) {
  double factorial_expr = 1.0;
  // factorial contribution of reactants
  for (int i = 0; i < current_reaction.reactant_types.size(); i++) {
    int nu_i = -1 * current_reaction.reactant_coefficients[i];
    int N_i0 = old_particle_numbers.at(current_reaction.reactant_types[i]);
    factorial_expr *= factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(
        N_i0, nu_i); // zeta = 1 (see @cite smith94c)
                     // since we only perform one reaction
                     // at one call of the function
  }
  // factorial contribution of products
  for (int i = 0; i < current_reaction.product_types.size(); i++) {
    int nu_i = current_reaction.product_coefficients[i];
    int N_i0 = old_particle_numbers.at(current_reaction.product_types[i]);
    factorial_expr *= factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(
        N_i0, nu_i); // zeta = 1 (see @cite smith94c)
                     // since we only perform one reaction
                     // at one call of the function
  }
  return factorial_expr;
}
} // namespace ReactionEnsemble
