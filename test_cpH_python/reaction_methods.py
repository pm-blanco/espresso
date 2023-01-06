class ReactionMethods:
    reactions=[]
    import random as rn
    import numpy as np
    particle_inside_exclusion_range_touched = False
    default_charges={}
    m_empty_p_ids_smaller_than_max_seen_particle=[]

    def __init__(self, **kwargs):
        self.kT=kwargs.pop("kT")
        self.exclusion_range=kwargs.pop("exclusion_range")
        self.seed=kwargs.pop("seed")
        self.exclusion_radius_per_type=kwargs.pop('exclusion_radius_per_type',{})
        self.np.random.seed(seed=self.seed)
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

        p_properties_old_state = self.make_reaction_attempt(current_reaction)
        if (self.particle_inside_exclusion_range_touched):
            # reject
            self.restore_system(system=system, p_properties_old_state=p_properties_old_state)
            return
        
        E_pot_new=system.analysis.energy()['non_bonded']
        bf = calculate_acceptance_probability(
            current_reaction, E_pot_old, E_pot_new, old_particle_numbers) # WIP
        
        exponential = self.np.exp(-1.0 / self.kT * (E_pot_new - E_pot_old))
        current_reaction.accumulator_potential_energy_difference_exponential.append(exponential)

        if (self.np.random.uniform(1) < bf):
            self.delete_hidden_particles(system=system, p_properties_old_state=p_properties_old_state)
            current_reaction.accepted_moves += 1
            E_pot_old = E_pot_new # Update the system energy
        else:
            # reject
            self.restore_system(system=system, p_properties_old_state=p_properties_old_state) 
        

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
            pid = self.create_particle(system=system, ptype=ptype) 
            self.check_exclusion_range(system=system, inserted_pid=pid) 
            p_properties_old_state['created_particle'].append({'pid':random_pid,
                                                              'type':ptype,
                                                              'charge':self.default_charges[ptype]})
        elif delta_n < 0:
          for _ in range(abs(delta_n)):
            random_pid=self.rn.choice(system.part.select(type=ptype))
            p_properties_old_state['hidden_particle'].append({'pid':random_pid,
                                                              'type':ptype,
                                                              'charge':self.default_charges[ptype]})
            self.check_exclusion_range(system=system, inserted_pid=pid)
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
            self.check_exclusion_range(system=system, inserted_pid=pid)
            self.hide_particle(system=system, pid=pid)
        else:
          # create additional product_types particles
          for _ in range(reaction.product_coefficients[index]):
            pid = self.create_particle(system=system, ptype=reaction.product_types[index])
            self.check_exclusion_range(system=system, inserted_pid=pid)
            p_properties_old_state['created_particle'].append({'pid':random_pid,
                                                              'type':ptype,
                                                              'charge':self.default_charges[ptype]})

      return p_properties_old_state
    
    def restore_system(self, **kwargs):

      system=kwargs.pop('system')
      p_properties_old_state=kwargs.pop('p_properties_old_state')

      # Restore the properties of changed_particles and hidden_particles

      for particle_info in p_properties_old_state['changed_particle']+p_properties_old_state['hidden_particle']:
        system.part.by_id(particle_info['pid']).type=particle_info['type']
        system.part.by_id(particle_info['pid']).q=particle_info['charge']

      # Destroy the created particles

      for particle_info in p_properties_old_state['created_particle']:
        system.by_id(particle_info['pid']).remove()

      return
    
    def hide_particle(self, **kwargs):
      system=kwargs.pop('system')
      pid=kwargs.pop('pid')
      system.part.by_id(pid).type=self.non_interacting_type
      system.part.by_id(pid).q=0
      return

    def delete_hidden_particles(self, **kwargs):
      system=kwargs.pop('system')
      for particle_info in kwargs.pop('p_properties_old_state'):
        system.part.by_id(particle_info['pid']).remove()
      return

                
    def create_particle(self, **kwargs):
      system=kwargs.pop('system')
      ptype=kwargs.pop('ptype')
      if len(self.m_empty_p_ids_smaller_than_max_seen_particle == 0):
        pid=system.part.highest_particle_id+1
      else:
        pid=min(self.m_empty_p_ids_smaller_than_max_seen_particle)
      
      # we use mass=1 for all particles, think about adapting this

      espresso_system.part.add(id=[pid], 
                              pos=[self.np.random.random((1, 3))[0] *self.np.copy(espresso_system.box_l)], 
                              type=[ptype], 
                              q=[self.default_charges[ptype]],
                              v=[self.np.random.normal(size=(1, 3))[0] *self.np.sqrt(self.kT)])

      return pid

    def check_exclusion_range(self, **kwargs):
      
      system=kwargs.pop('system')
      inserted_particle=system.part.by_id(kwargs.pop('inserted_pid'))
      if self.exclusion_radius_per_type:
        if self.exclusion_radius_per_type[inserted_particle.type] == 0:
          return

      # NOTE: Missing the feature of searching short range neighbours
  
      for particle in system.part.all():
        
        if particle.id == inserted_pid:
          continue

        if (self.exclusion_radius_per_type[inserted_particle.type] == 0)  or  (self.exclusion_radius_per_type[particle.type] == 0):
          excluded_distance = self.exclusion_range
        else:
          excluded_distance = self.exclusion_radius_per_type[inserted_particle.type] +  self.exclusion_radius_per_type[particle.type]

        d_min=system.distance(particle, inserted_particle)
        if d_min < excluded_distance:
          self.particle_inside_exclusion_range_touched=True
          return

      return


"""


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
        self.accumulator_potential_energy_difference_exponential=[]
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

            