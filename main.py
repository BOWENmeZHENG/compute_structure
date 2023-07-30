
from cProfile import label
from operator import index
from processing import *
import numpy as np
import matplotlib.pyplot as plt

"""
Compute bond length and tilting angle of MD simulations
"""

BOND_LENGTH_MEAN, BOND_lENGTH_STD = 2.23, 0.114
ANGLE_MEAN, ANGLE_STD = 118.62, 10.85

trj = parse_trj("torch_3_4A_20k_5.dump")

bond_length = []
angle_Mg_O_C = []

for frame in trj:
    bl = comp_ave_bond_len(frame)
    bond_length.append(bl)
    angle = comp_ave_tilt_angle(frame)
    angle_Mg_O_C.append(angle)

time = np.arange(0, len(bond_length), 1) * 100 * 0.0005

fig, ax1= plt.subplots(1, 1, figsize=(8, 6))
ax1.set_title("SNAP potential trained with 600 K data", fontsize=18)
ax1.plot(time, bond_length, linewidth=1, label="SNAP")
ax1.axhline(BOND_LENGTH_MEAN, color="red", linewidth=1, label="DFT")
ax1.axhspan(BOND_LENGTH_MEAN - BOND_lENGTH_STD, BOND_LENGTH_MEAN + BOND_lENGTH_STD, alpha=0.3, color='red')
ax1.set_xlabel("time (ps)", fontsize=18)
ax1.set_ylabel(r"Mg-$\rm O_{CO2}$ bond length ($\rm \AA$)", fontsize=18)
ax1.tick_params(labelsize=16)
ax1.legend(fontsize=18)

plt.show()

# indices = [i for i, bl in enumerate(bond_length) if bl >3]
# long_bonds = [bond_length[index] for index in indices]
# angles = [angle_Mg_O_C[index] for index in indices]
# print(len(indices))
# print(angles)
# plt.figure(figsize=(8, 6))
# # plt.scatter(long_bonds[1:], angles[1:])
# plt.scatter(indices[1:], angles[1:])
# plt.show()

fig, ax1= plt.subplots(1, 1, figsize=(8, 6))
ax1.set_title("SNAP potential trained with 600 K data", fontsize=18)
ax1.plot(time, angle_Mg_O_C, linewidth=1, label="SNAP")
ax1.axhline(ANGLE_MEAN, color="red", linewidth=1, label="DFT")
ax1.axhspan(ANGLE_MEAN - ANGLE_STD, ANGLE_MEAN + ANGLE_STD, alpha=0.3, color='red')
ax1.set_ylim(bottom=0, top=180)
ax1.set_xlabel("time (ps)", fontsize=18)
ax1.set_ylabel(r"tilting angle $\rm Mg-O_{CO2}-C_{CO2}$ ($\rm \degree $)", fontsize=18)
ax1.tick_params(labelsize=16)
ax1.legend(fontsize=18)

plt.show()
