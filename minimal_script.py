import espressomd
import espressomd.accumulators
import espressomd.observables
import numpy as np

NPs_number=10
system = espressomd.System(box_l=[10]*3)
system.time_step = 0.1
system.cell_system.skin = 0.4

for _ in range(int(NPs_number)):  
    system.part.add(pos  = np.random.random((1, 3))[0]*10)

obs = espressomd.observables.ContactTimes(ids=system.part.all().id)
accumulator = espressomd.accumulators.TimeSeries(obs=obs, 
                                                delta_N=1)