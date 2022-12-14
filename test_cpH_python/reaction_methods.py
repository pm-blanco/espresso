class ReactionMethods:
    reactions=[]
    import random as rn
    import numpy as np
    particle_inside_exclusion_range_touched = False
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
            self.restore_system()
            return
        
        E_pot_new=system.analysis.energy()['non_bonded']
        bf = calculate_acceptance_probability(
            current_reaction, E_pot_old, E_pot_new, old_particle_numbers)
        
        #exponential = self.np.exp(-1.0 / self.kT * (E_pot_new - E_pot_old))
        #current_reaction.accumulator_potential_energy_difference_exponential(
        #    exponential);

        if (self.np.random.uniform(1) < bf):
            change_tracker.delete_hidden_particles()
            current_reaction.accepted_moves += 1
            E_pot_old = E_pot_new # Update the system energy
        else:
            # reject
            change_tracker.restore_original_state()
        

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

"""
ReactionAlgorithm::make_reaction_attempt(SingleReaction const &reaction) {
  // create or hide particles of types with corresponding types in reaction
  auto const n_product_types = reaction.product_types.size();
  auto const n_reactant_types = reaction.reactant_types.size();
  auto const get_random_p_id_of_type = [this](int type) {
    auto const random_index = i_random(number_of_particles_with_type(type));
    return get_random_p_id(type, random_index);
  };
  ParticleChangeRecorder tracker{[this](int p_id) { delete_particle(p_id); }};
  for (int i = 0; i < std::min(n_product_types, n_reactant_types); i++) {
    auto const n_product_coef = reaction.product_coefficients[i];
    auto const n_reactant_coef = reaction.reactant_coefficients[i];
    // change std::min(reactant_coefficients(i),product_coefficients(i)) many
    // particles of reactant_types(i) to product_types(i)
    auto const type = reaction.reactant_types[i];
    for (int j = 0; j < std::min(n_product_coef, n_reactant_coef); j++) {
      auto const p_id = get_random_p_id_of_type(type);
      tracker.save_changed_particle({p_id, type, charges_of_types[type]});
      replace_particle(p_id, reaction.product_types[i]);
    }
    // create product_coefficients(i)-reactant_coefficients(i) many product
    // particles iff product_coefficients(i)-reactant_coefficients(i)>0,
    // iff product_coefficients(i)-reactant_coefficients(i)<0, hide this number
    // of reactant particles
    auto const delta_n = n_product_coef - n_reactant_coef;
    if (delta_n > 0) {
      auto const type = reaction.product_types[i];
      for (int j = 0; j < delta_n; j++) {
        auto const p_id = create_particle(type);
        check_exclusion_range(p_id);
        tracker.save_created_particle(p_id);
      }
    } else if (delta_n < 0) {
      auto const type = reaction.reactant_types[i];
      for (int j = 0; j < -delta_n; j++) {
        auto const p_id = get_random_p_id_of_type(type);
        tracker.save_hidden_particle({p_id, type, charges_of_types[type]});
        check_exclusion_range(p_id);
        hide_particle(p_id);
      }
    }
  }
  // create or hide particles of types with noncorresponding replacement types
  for (auto i = std::min(n_product_types, n_reactant_types);
       i < std::max(n_product_types, n_reactant_types); i++) {
    if (n_product_types < n_reactant_types) {
      auto const type = reaction.reactant_types[i];
      // hide superfluous reactant_types particles
      for (int j = 0; j < reaction.reactant_coefficients[i]; j++) {
        auto const p_id = get_random_p_id_of_type(type);
        tracker.save_hidden_particle({p_id, type, charges_of_types[type]});
        check_exclusion_range(p_id);
        hide_particle(p_id);
      }
    } else {
      // create additional product_types particles
      for (int j = 0; j < reaction.product_coefficients[i]; j++) {
        auto const p_id = create_particle(reaction.product_types[i]);
        check_exclusion_range(p_id);
        tracker.save_created_particle(p_id);
      }
    }
  }

  return tracker;
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

            