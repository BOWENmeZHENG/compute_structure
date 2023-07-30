import numpy as np
import utils

ATOM = 327
USELESS_LINE = 8


def parse_trj(trj_name: str) -> list:
    '''
    Parse the entire trajectory, store frames in a list
    '''
    trj = []
    with open(trj_name) as f_trj:
        line = f_trj.readline()
        while True:
            frame = []
            for i in range(USELESS_LINE):
                line = f_trj.readline()
            for i in range(ATOM):
                line = f_trj.readline()
                line_list = [float(num_str) for num_str in line.split()]
                line_list[0] = int(line_list[0])
                line_list[1] = int(line_list[1])
                frame.append(line_list)
            line = f_trj.readline()
            trj.append(frame)
            if line == "":
                break
    return trj



def bonded_O_Mg(O_CO2_list, C_CO2, Mg_list):
    O_bonded = None
    Mg_bonded = None
    bond_length = None
    dist_1, dist_2 = [], []
    first_O = np.array(O_CO2_list[0][2:5])
    second_O = np.array(O_CO2_list[1][2:5])
    C_xyz = np.array(C_CO2[2:5])
    first_O, C_xyz, second_O = utils.checkwrap(first_O, C_xyz, second_O)

    for Mg in Mg_list:
        dist_1.append(utils.comp_dist(first_O, np.array(Mg[2:5])))
        dist_2.append(utils.comp_dist(second_O, np.array(Mg[2:5])))
    min_dist_1 = min(dist_1)
    min_dist_2 = min(dist_2)
    if min_dist_1 < min_dist_2:
        O_bonded = first_O
        Mg_index = dist_1.index(min_dist_1)
        Mg_bonded = Mg_list[Mg_index]
        bond_length = min_dist_1
    else:
        O_bonded = second_O
        Mg_index = dist_2.index(min_dist_2)
        Mg_bonded = Mg_list[Mg_index]
        bond_length = min_dist_2

    return O_bonded, np.array(Mg_bonded[2:5]), bond_length



def comp_ave_bond_len(frame: list):
    O_CO2_list, C_CO2, Mg_list = utils.list_atoms(frame)
    *_, bond_length = bonded_O_Mg(O_CO2_list, C_CO2, Mg_list)
    return bond_length

def comp_ave_tilt_angle(frame: list):
    O_CO2_list, C_CO2, Mg_list = utils.list_atoms(frame)
    O_bonded, Mg_bonded, _ = bonded_O_Mg(O_CO2_list, C_CO2, Mg_list)
    tilt_angles = utils.comp_angle(Mg_bonded, O_bonded, np.array(C_CO2[2:5]))
    return tilt_angles