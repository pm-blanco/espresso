import espressomd
import espressomd.accumulators
import espressomd.observables
import numpy as np
import json 


def calculate_contact_times(accumulator_ids1, accumulator_ids2, contact_threshold, system):
    """
    Calculates the contact times from the time series of the coordinates stored in `accumulator_ids1`  and `accumulator_ids2`.

    Args:
        accumulator_ids1(`espressomd.accumulators.TimeSeries`): accumulator with the time series of the first set of particles.
        accumulator_ids2(`dict`): accumulator with the time series of the second set of particles.
        contact_threshold(`float`): threshold distance under which particles are considered to be in contact.
        system(`espressomd.System`): instance of the espresso simulation system.
    """
    # Get the coordinates 
    coords_ids = accumulator_ids1.time_series()
    coords_target = accumulator_ids2.time_series()
    # Check for the initial contacts (if any)
    first_coord_set_ids = coords_ids[0]
    first_coord_set_target = coords_target[0]
    n_coords = len(coords_ids)
    n_ids1 = len(first_coord_set_ids)
    n_ids2 = len(first_coord_set_target)
    contacts = np.zeros((n_ids1, n_ids2))
    time_first_contact = np.zeros((n_ids1, n_ids2))
    # Calculate contact times:
    frame_time = 0
    contact_times = []
    for n_coord in range(n_coords):
        coord_set_ids = coords_ids[n_coord]
        coord_set_target = coords_target[n_coord]
        for index1 in range(n_ids1):
            for index2 in range(n_ids2):
                if index1 == index2:
                    continue
                dist = calculate_minimum_image_distance(coord_set_ids[index1],
                                                        coord_set_target[index2], 
                                                        system.box_l)
                if dist < contact_threshold:   # They are in contact now
                    # but they were not contact before!
                    if not contacts[index1][index2]:
                        contacts[index1][index2] = True
                        time_first_contact[index1][index2] = frame_time
                else:  # They are not in contact
                    if contacts[index1][index2]:  # but they were in contact before! 
                        # Calculate the total contact time
                        contact_time = frame_time - \
                            time_first_contact[index1][index2] - \
                            system.time_step
                        if (contact_time < system.time_step):  # Avoid rounding errors
                            contact_time = 0
                        contacts[index1][index2] = False
                        contact_times.append(float(contact_time))
        frame_time += system.time_step
    return contact_times


def calculate_minimum_image_distance(coords1, coords2, size_box):
    """
    Calculates the minimum distance between two particles in the simulation box,
    taking into account the minimum image convention and periodic boundary conditions.

    Args:
        coords1(`list`): array with the cartesian coordinates of the first particle.
        coords2(`list`): array with the cartesian coordinates of the second particle.
        size_box(`list`): xyz size of the simulation box.

    Returns:
        dist(`float`): minimum distance between the particles
    """
    import numpy as np
    # Sanity check, check that all input variables have the same size
    if (len(coords1) != len(coords2)) or (len(coords2) != len(size_box)):
        raise ValueError(
            f"All input arrays must have the same size. The sizes are: coords1 = {len(coords1)}, coords2 = {len(coords2)} and size_box = {len(size_box)}")

    n_axes = len(coords1)
    coords = {"particle1": coords1,
              "particle2": coords2}

    # Fold particles
    for particle in coords.keys():
        for axis in range(n_axes):
            coords[particle][axis] = coords[particle][axis] % size_box[axis]

    # Calculate minimum distance
    dist2 = 0
    for axis in range(n_axes):
        dist = abs(coords["particle1"][axis] - coords["particle2"][axis])
        if dist > size_box[axis] * 0.5:
            dist -= size_box[axis]
        dist2 += dist**2
    return np.sqrt(dist2)


# Switch to use the Python or the C++ implementation
mode = "c++"  # supported modes are 'python' and 'c++'

# System: LJ fluid
cutoff = 3
offset = 1
epsilon = 1
sigma = 1
contact_threshold = cutoff + offset
N = 50
seed = 42
volume_fraction = 0.01
effective_radi = (sigma + offset) / 2
# Simulation parameters
N_steps = 100000
box_l = (4. / 3. * N * np.pi * effective_radi**3 / volume_fraction)**(1. / 3.)
system = espressomd.System(box_l=[box_l] * 3)
system.time_step = 0.01
system.cell_system.skin = 1
ids = range(N)
rng = np.random.default_rng(seed)
system.part.add(pos=rng.random((N, 3)) * np.copy(system.box_l))

# relax the system
system.integrator.set_steepest_descent(f_max=1e-3, 
                                       gamma=1, 
                                       max_displacement=0.1)
system.integrator.run(1000)

# setup LD integration
system.thermostat.set_langevin(kT=1, 
                               gamma=1, 
                               seed=seed)
system.integrator.set_vv() 

# Setup LJ interaction
system.non_bonded_inter[0, 0].lennard_jones.set_params(epsilon=epsilon,
                                                       sigma=sigma,
                                                       cutoff=cutoff,
                                                       offset=offset,
                                                       shift="auto")


# Setup the observables to track the contact time
if mode == "c++":
    obs = espressomd.observables.ContactTimes(ids=ids,
                                              target_ids=ids,
                                              contact_threshold=contact_threshold)

elif mode == "python":
    # in python, the contact time needs to be calculated a posteriori from the time series of the trajectories
    obs = espressomd.observables.ParticlePositions(ids=ids)
    obs_target_ids = espressomd.observables.ParticlePositions(ids=ids)
    accumulator_target_ids = espressomd.accumulators.TimeSeries(obs=obs_target_ids, 
                                                                delta_N=1)
    system.auto_update_accumulators.add(accumulator_target_ids)    

accumulator = espressomd.accumulators.TimeSeries(obs=obs, 
                                                 delta_N=1)
system.auto_update_accumulators.add(accumulator)

# Run the simulation
system.integrator.run(N_steps)

# Get the contact time
if mode == "c++":
    contact_times = list(obs.contact_times_series())
elif mode == "python":
    contact_times = calculate_contact_times(accumulator_ids1=accumulator,
                                            accumulator_ids2=accumulator_target_ids,
                                            contact_threshold=contact_threshold,
                                            system=system)


# Output the data as a json dictionary
with open(f"contact_times_{mode}.json", "w") as outfile: 
    json.dump({"contact_times": contact_times}, outfile)
