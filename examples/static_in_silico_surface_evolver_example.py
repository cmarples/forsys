import sys
sys.path.append('..')
sys.path.insert(0, '../../forsys')
import forsys as fs
import os

DATA_FOLDER = os.path.join("data", "in_silico")
RESULTS_FOLDER = os.path.join("results")

# This is only necessary if you wish to create outputs
# if not os.path.exists(RESULTS_FOLDER):
#     os.makedirs(RESULTS_FOLDER)

lattice = fs.surface_evolver.SurfaceEvolver(os.path.join(DATA_FOLDER,
                                                           f"step_{24}.dmp"))

frames = {}
frames[0] = fs.frames.Frame(0,
                            lattice.vertices,
                            lattice.edges, 
                            lattice.cells,
                            gt=True)
forsys = fs.ForSys(frames)

forsys.build_force_matrix(when=0)
forsys.solve_stress(when=0, allow_negatives=False)

forsys.build_pressure_matrix(when=0)
forsys.solve_pressure(when=0, method="lagrange_pressure")

fig, ax = fs.plot.plot_inference(forsys.frames[0],
                                 ground_truth=True,
                                 normalized="max",
                                 mirror_y=False,
                                 colorbar=False,
                                 pressure=False)

fig, ax = fs.plot.plot_inference(forsys.frames[0],
                                 normalized="max",
                                 mirror_y=False,
                                 colorbar=False,
                                 pressure=False)

fig, ax = fs.plot.plot_inference(forsys.frames[0],
                                 normalized="max",
                                 mirror_y=False,
                                 colorbar=False,
                                 pressure=True)

