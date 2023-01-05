class ReactionMethods:
    reactions=[]
    import random as rn
    import numpy as np
    particle_inside_exclusion_range_touched = False
    default_charges={}
    def __init__(self, **kwargs):
        self.kT=kwargs.pop("kT")
        self.exclusion_range=kwargs.pop("exclusion_range")
        self.seed=kwargs.pop("seed")
        ReactionMethods.np.random.seed(seed=self.seed)
        return

    def add_reaction(self, **kwargs):
        reaction=SingleReaction(reactant_types=kwargs.pop('reactant_types'),
                                reactant_coefficients=kwargs.pop('reactant_coefficients'),
                                product_types=kwargs.pop('product_types'),
                                product_coefficients=kwargs.pop('product_coefficients'),
                                Gamma=kwargs.pop('Gamma'),)
        for ptype in kwargs.pop('default_charges').keys():
          self.default_charges[ptype]=default_charges[ptype]

        self.reactions.append(reaction)
        return
    
    def do_reaction(self, **kwargs):
        system=kwargs.pop('system')
        E_pot_old=system.analysis.energy()['non_bonded']
        for _ in range(kwargs.pop('reaction_steps')):
            reaction=self.rn.choice(self.reactions)
            self.generic_oneway_reaction(reaction=reaction,
                                         E_pot_old=E_pot_old,
                                        system=system)

    def generic_oneway_reaction(self,**kwargs):

        current_reaction=kwargs.pop('reaction')
        system=kwargs.pop('system')
        E_pot_old=kwargs.pop('E_pot_old')
        current_reaction.trial_moves+=1
        if not self.all_reactant_particles_exist(reaction=current_reaction,
                                                    system=system):
            return

        old_particle_numbers = self.save_old_particle_numbers(reaction=current_reaction,
                                                                system=system)

        change_tracker = self.make_reaction_attempt(current_reaction) # WIP
        if (self.particle_inside_exclusion_range_touched):
            self.restore_system() # TODO
            return
        
        E_pot_new=system.analysis.energy()['non_bonded']
        bf = calculate_acceptance_probability(
            current_reaction, E_pot_old, E_pot_new, old_particle_numbers) # WIP
        
        #exponential = self.np.exp(-1.0 / self.kT * (E_pot_new - E_pot_old))
        #current_reaction.accumulator_potential_energy_difference_exponential(
        #    exponential);

        if (self.np.random.uniform(1) < bf):
            change_tracker.delete_hidden_particles()
            current_reaction.accepted_moves += 1
            E_pot_old = E_pot_new # Update the system energy
        else:
            # reject
            self.restore_system() # TODO
        

        return

    def all_reactant_particles_exist(self,**kwargs):
        reaction=kwargs.pop('reaction')
        system=kwargs.pop('system')
        for r_type in reaction.reactant_types:
            index=reaction.reactant_types.index(r_type)
            if (system.number_of_particles(type=r_type) < reaction.reactant_coefficients[index]):
                return False

        return True

    def save_old_particle_numbers(self, **kwargs):
        old_particle_numbers={}
        reaction=kwargs.pop('reaction')
        system=kwargs.pop('system')
        for r_type in reaction.reactant_types+reaction.product_types:
            old_particle_numbers[r_type]=system.number_of_particles(type=r_type)
        return old_particle_numbers
  
    def make_reaction_attempt(self, **kwargs): 
  
      reaction = kwargs.pop('reaction')

      # create or hide particles of types with corresponding types in reaction
    
      minimum_number_of_types=min(len(reaction.reactant_types), len(reaction.product_types))
      p_properties_old_state={'changed_particle': [],
                              'created_particle': [],
                              'hidden_particle': []}

      for index  in range(minimum_number_of_types):

        # change minimum_number_of_types many particles 
        # of reactant_types(i) to product_types(i)
        
        min_stechiometic_coeff=min(reaction.reactant_coefficients[index], 
                                    reaction.product_coefficients[index])
        ptype=reaction.reactant_types[index]
        for _ in range(min_stechiometic_coeff):
          random_pid=self.rn.choice(system.part.select(type=ptype))
          p_properties_old_state['changed_particle'].append({'pid':random_pid,
                                                              'type':ptype,
                                                              'charge':self.default_charges[ptype]})

        # create product_coefficients(i)-reactant_coefficients(i) many product
        # particles iff product_coefficients(i)-reactant_coefficients(i)>0,
        # iff product_coefficients(i)-reactant_coefficients(i)<0, hide this number
        # of reactant particles

        delta_n=reaction.product_coefficients[index]-reaction.reactant_coefficients[index]

        if delta_n > 0:
          ptype=reaction.product_types[index]
          for _ in range(delta_n):
            pid = self.create_particle(p_type) # WIP
            self.check_exclusion_range(pid) # WIP
            p_properties_old_state['created_particle'].append({'pid':random_pid,
                                                              'type':ptype,
                                                              'charge':self.default_charges[ptype]})
        elif delta_n < 0:
          for _ in range(abs(delta_n)):
            random_pid=self.rn.choice(system.part.select(type=ptype))
            p_properties_old_state['hidden_particle'].append({'pid':random_pid,
                                                              'type':ptype,
                                                              'charge':self.default_charges[ptype]})
            self.check_exclusion_range(p_id); # WIP
            self.hide_particle(system=system, pid=pid)

    # create or hide particles of types with noncorresponding replacement types
      maximum_number_of_types=max(len(reaction.reactant_types), len(reaction.product_types))
      for index in range(minimum_number_of_types,maximum_number_of_types):
        if len(reaction.product_types) < len(reaction.reactant_types):
          # hide superfluous reactant_types particles
          ptype=reaction.reactant_types[index]
          for _ in range(reaction.reactant_coefficients[index]):
            pid = self.rn.choice(system.part.select(type=ptype))
            p_properties_old_state['hidden_particle'].append({'pid':random_pid,
                                                              'type':ptype,
                                                              'charge':self.default_charges[ptype]})
            self.check_exclusion_range(pid) # WIP
            self.hide_particle(system=system, pid=pid)
        else:
          # create additional product_types particles
          for _ in range(reaction.product_coefficients[index]):
            pid = self.create_particle(reaction.product_types[index]) # WIP
            self.check_exclusion_range(pid) # WIP
            p_properties_old_state['created_particle'].append({'pid':random_pid,
                                                              'type':ptype,
                                                              'charge':self.default_charges[ptype]})


    def hide_particle(self, **kwargs):
      system=kwargs.pop('system')
      pid=kwargs.pop('pid')
      system.by_id(pid).type=self.non_interacting_type
      system.by_id(pid).q=0
      return
"""


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

  // we use mass=1 for all particles, think about adapting this
  auto const new_pos = get_random_position_in_box();
  mpi_make_new_particle(p_id, new_pos);
  move_particle(p_id, new_pos, std::sqrt(kT));
  set_particle_type(p_id, desired_type);
#ifdef ELECTROSTATICS
  set_particle_q(p_id, charges_of_types[desired_type]);
#endif
  return p_id;
}

void ReactionAlgorithm::move_particle(int p_id, Utils::Vector3d const &new_pos,
                                      double velocity_prefactor) {
  mpi_set_particle_pos(p_id, new_pos);
  // create random velocity vector according to Maxwell-Boltzmann distribution
  Utils::Vector3d vel;
  vel[0] = velocity_prefactor * m_normal_distribution(m_generator);
  vel[1] = velocity_prefactor * m_normal_distribution(m_generator);
  vel[2] = velocity_prefactor * m_normal_distribution(m_generator);
  set_particle_v(p_id, vel);
}

/**
 * Check if the inserted particle is too close to neighboring particles.
 */
void ReactionAlgorithm::check_exclusion_range(int inserted_particle_id) {

  auto const &inserted_particle = get_particle_data(inserted_particle_id);

  /* Check the exclusion radius of the inserted particle */
  if (exclusion_radius_per_type.count(inserted_particle.type()) != 0) {
    if (exclusion_radius_per_type[inserted_particle.type()] == 0.) {
      return;
    }
  }

  std::vector<int> particle_ids;
  if (neighbor_search_order_n) {
    particle_ids = get_particle_ids();
    /* remove the inserted particle id */
    particle_ids.erase(std::remove(particle_ids.begin(), particle_ids.end(),
                                   inserted_particle_id),
                       particle_ids.end());
  } else {
    particle_ids = mpi_get_short_range_neighbors(inserted_particle.id(),
                                                 m_max_exclusion_range);
  }

  /* Check if the inserted particle within the exclusion radius of any other
   * particle */
  for (auto const &particle_id : particle_ids) {
    auto const &p = get_particle_data(particle_id);
    double excluded_distance;
    if (exclusion_radius_per_type.count(inserted_particle.type()) == 0 ||
        exclusion_radius_per_type.count(p.type()) == 0) {
      excluded_distance = exclusion_range;
    } else if (exclusion_radius_per_type[p.type()] == 0.) {
      continue;
    } else {
      excluded_distance = exclusion_radius_per_type[inserted_particle.type()] +
                          exclusion_radius_per_type[p.type()];
    }

    auto const d_min =
        box_geo.get_mi_vector(p.pos(), inserted_particle.pos()).norm();

    if (d_min < excluded_distance) {
      particle_inside_exclusion_range_touched = true;
      break;
    }
  }
}

"""


class SingleReaction:
    def __init__(self,**kwargs):
        self.reactant_types=kwargs.pop('reactant_types')
        self.reactant_coefficients=kwargs.pop('reactant_coefficients')
        self.product_types=kwargs.pop('product_types')
        self.product_coefficients=kwargs.pop('product_coefficients')
        self.Gamma=kwargs.pop('Gamma')
        self.accepted_moves=0
        self.trial_moves=0
        return

class ConstantpHEnsemble(ReactionMethods):
    def __init__(self, **kwargs):
        self=ReactionMethods(kT=kwargs.pop("kT"),
                              exclusion_range=kwargs.pop("exclusion_range"),
                             seed=kwargs.pop("seed"))
        self.constant_pH=kwargs.pop("constant_pH")
        return

    def add_acidbase_reaction(self, **kwargs):
        Gamma=10**-(kwargs.pop('pKa'))
        # Add reaction in the forwards direction
        ReactionMethods().add_reaction(reactant_types=kwargs.pop('reactant_types'),
                                    reactant_coefficients=[1],
                                    product_types=kwargs.pop('product_types'),
                                    product_coefficients=[1,1],
                                    Gamma=Gamma
                                    )
        # Add reaction in the reverse direction
        ReactionMethods().add_reaction(reactant_types=kwargs.pop('product_types'),
                                    reactant_coefficients=[1,1],
                                    product_types=kwargs.pop('reactant_types'),
                                    product_coefficients=[1],
                                    Gamma=Gamma
                                    )

            