import numpy as np

C = 13.992


def checkwrap(xyz_a, xyz_b, xyz_c, pbc_th=6):
    if abs(xyz_a[2] - xyz_b[2]) > pbc_th \
    or abs(xyz_b[2] - xyz_c[2]) > pbc_th \
    or abs(xyz_c[2] - xyz_a[2]) > pbc_th:
        if xyz_a[2] < C / 2 and xyz_b[2] > C / 2 and xyz_c[2] > C / 2:
            xyz_a[2] += C
        if xyz_a[2] > C / 2 and xyz_b[2] < C / 2 and xyz_c[2] < C / 2:
            xyz_a[2] -= C
        if xyz_b[2] < C / 2 and xyz_c[2] > C / 2 and xyz_a[2] > C / 2:
            xyz_b[2] += C
        if xyz_b[2] > C / 2 and xyz_c[2] < C / 2 and xyz_a[2] < C / 2:
            xyz_b[2] -= C
        if xyz_c[2] < C / 2 and xyz_b[2] > C / 2 and xyz_c[2] > C / 2:
            xyz_c[2] += C
        if xyz_c[2] > C / 2 and xyz_b[2] < C / 2 and xyz_c[2] < C / 2:
            xyz_c[2] -= C
    return xyz_a, xyz_b, xyz_c

def comp_dist(xyz_1, xyz_2):
    return np.linalg.norm(xyz_2 - xyz_1)


def comp_angle(xyz_a, xyz_b, xyz_c):
    xyz_a, xyz_b, xyz_c = checkwrap(xyz_a, xyz_b, xyz_c)
    ba = xyz_a - xyz_b
    bc = xyz_c - xyz_b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)


def list_atoms(frame: list):
    O_CO2_list = []  # list of CO2 O atoms
    C_CO2 = None
    Mg_list = []  # list of Mg atoms
    for atom in frame:
        if atom[1] == 6:
            O_CO2_list.append(atom)
        if atom[1] == 5:
            C_CO2 = atom
        if atom[1] == 1:
            Mg_list.append(atom)
    return O_CO2_list, C_CO2, Mg_list 

def sort_lists_by_index(atom_list):
    return sorted(atom_list, key=lambda x: x[0])
    
