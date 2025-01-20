import espressomd
import espressomd.accumulators
import espressomd.observables
import numpy as np

NPs_number=2
system = espressomd.System(box_l=[10]*3)
system.time_step = 0.01
system.cell_system.skin = 1


system.part.add(pos  = [0,0,0])
system.part.add(pos  = [2,2,2])

system.part.add(pos  = [1,1,1])
system.part.add(pos  = [6,6,6])

obs = espressomd.observables.ContactTimes(ids=[0,1],
                                          target_ids=[2,3],
                                          contact_threshold=4)
accumulator = espressomd.accumulators.TimeSeries(obs=obs, 
                                                delta_N=1)
system.thermostat.set_langevin(kT=1, 
                                gamma=1, 
                                seed=1)

system.auto_update_accumulators.add(accumulator)

system.non_bonded_inter[0,0].lennard_jones.set_params(epsilon = 1,
                                                    sigma   = 1,
                                                    cutoff  = 3,
                                                    offset  = 1,
                                                    shift   = "auto")


system.integrator.run(1000)
res=obs.contact_times_series()
current_times = obs.get_instantaneous_contact_times()
print(type(res))
# Get indices where res is non-zero
indices, = np.nonzero(res)  # np.nonzero returns a tuple
print(type(indices))
# Access elements at those indices
elements = res[indices]
print(elements)
obs.clean_contact_times()
res=obs.contact_times_series()
print(res)
current_times = obs.get_instantaneous_contact_times()
print(current_times)
