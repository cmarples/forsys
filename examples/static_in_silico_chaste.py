# Script for using Chaste output as ForSys input

import sys
sys.path.append('..')
import os
import forsys as fs

CHASTE_OUTPUT_NAME = "X_Furrow_K_0p1"
FRAME_FILE = os.path.join("results_from_time_0", "results_5000.vtu")

CHASTE_OUTPUT_FILE = os.path.join("in_silico_chaste", CHASTE_OUTPUT_NAME)
CHASTE_OUTPUT_FOLDER = os.path.join("data", CHASTE_OUTPUT_FILE)
FRAME_FILE_PATH = os.path.join(CHASTE_OUTPUT_FOLDER, FRAME_FILE)

# Method used to calculate versors
REPRESENTATION_METHOD = "straight"

# Load Chaste frame into a lattice
lattice = fs.chaste.Chaste(CHASTE_OUTPUT_FOLDER, FRAME_FILE_PATH)

# Assign frames and create main ForSys object
frames = {}
frames[0] = fs.frames.Frame(0,
                            lattice.vertices,
                            lattice.edges, 
                            lattice.cells,
                            gt=True)
forsys = fs.ForSys(frames)

# Build and solve the system of equations for the force
forsys.build_force_matrix(when=0, representation_method=REPRESENTATION_METHOD)
forsys.solve_stress(when=0, allow_negatives=False)

# Build and solve the system of equations for the pressure
forsys.build_pressure_matrix(when=0)
forsys.solve_pressure(when=0, method="lagrange_pressure")

max_force = None

# Plot with ground truth
fig1, ax1 = fs.plot.plot_inference(forsys.frames[0],
                                   ground_truth=True,
                                   normalized="max",
                                   maxForce=max_force,
                                   mirror_y=False,
                                   colorbar=True,
                                   pressure=False)

# Plot with inferred
fig2, ax2 = fs.plot.plot_inference(forsys.frames[0],
                                 normalized="max",
                                 maxForce=max_force,
                                 mirror_y=False,
                                 colorbar=True,
                                 pressure=False)

import matplotlib.pyplot as plt
plt.figure(figsize=(10, 10), dpi=600)
fs.plot.plot_mesh(forsys.frames[0], ax2, plot_tjs=True, plot_versors=True, representation_method=REPRESENTATION_METHOD)






