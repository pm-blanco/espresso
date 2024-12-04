import espressomd
import espressomd.accumulators
import espressomd.observables
import numpy as np

NPs_number=2
system = espressomd.System(box_l=[10]*3)
system.time_step = 0.001
system.cell_system.skin = 1


system.part.add(pos  = [0,0,0])
system.part.add(pos  = [2,2,2])

system.part.add(pos  = [5,5,5])
system.part.add(pos  = [6,6,6])

obs = espressomd.observables.ContactTimes(ids=[0,1],
                                          max_z=2)
accumulator = espressomd.accumulators.TimeSeries(obs=obs, 
                                                delta_N=1)
system.thermostat.set_langevin(kT=1, 
                                gamma=1, 
                                seed=1)
system.non_bonded_inter[0,0].lennard_jones.set_params(epsilon = 1,
                                                    sigma   = 1,
                                                    cutoff  = 2**(1./6.),
                                                    offset  = 1,
                                                    shift   = "auto",
                                                    )

system.integrator.run(1)
obs.calculate()

system.integrator.run(10000)
obs.calculate()