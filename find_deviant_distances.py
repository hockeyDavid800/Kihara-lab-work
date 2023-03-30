import numpy
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.vectors import Vector, calc_dihedral, calc_angle
import matplotlib.pyplot as plt
import numpy as np
import pathlib
from collections import defaultdict


# this program will run on each RNA
rna_angles_atoms = [
        ["O3'", "P", "O5'", "C5'"],  # alpha
        ["P", "O5'", "C5'", "C4'"],  # beta
        ["O5'", "C5'", "C4'", "C3'"],  # gamma
        ["C5'", "C4'", "C3'", "O3'"],  # delta
        ["C4'", "C3'", "O3'", "P"],  # epsilon
        ["C3'", "O3'", "P", "O5'"],  # zeta
        ["O4'", "C1'", "N1", "C2"],  # pyrimidine
        ["O4'", "C1'", "N9", "C4"]  # purine
    ]

hinge_atoms_base = [
        ["OP1", "P", "OP2"],
        ["P", "O5'", "C5'"],
        ["O5'", "C5'", "C4'"],
        ["C5'", "C4'", "O4'"],
        ["C5'", "C4'", "C3'"],
        ["O4'", "C4'", "C3'"],
        ["C4'", "C3'", "C2'"],
        ["O3'", "C3'", "C2'"],
        ["C3'", "C2'", "O2'"],
        ["C3'", "C2'", "C1'"],
        ["C2'", "C1'", "O4'"],
        ["C4'", "C3'", "O3'"],
        ["O3'", "C3'", "C2'"],
        ["C3'", "C2'", "O2'"],
        ["O2'", "C2'", "C1'"],
        ["C1'", "O4'", "C4'"]
                            ]
#specials below:
hinge_atoms_pu = [
                         ["C1'", "N9", "C8"],
                        ["N9", "C8", "N7"],
                        ["C8", "N7", "C5"],
                        ["N7", "C5", "C6"],
                        ["C5", "C6", "N1"],
                        ["C6", "N1", "C2"],
                        ["N1", "C2", "N3"],
                        ["C2", "N3", "C4"],
                        ["N3", "C4" , "C5"],
                        ["N9", "C4", "C5"],
                        ["N1", "C2", "N2'"],
                        ["N2", "C2", "N3"],
                        ["C5", "C6", "O6"],
                        ["O6", "C6", "N1"],
                        ["C4", "C5", "C6"],
                        ["O4'", "C1'", "N9"]
                                                ]
hinge_atoms_py = [
                       ["C1'", "N1", "C2"],
                       ["N1", "C2", "C3"],
                       ["C2", "C3", "C4"],
                       ["C3", "C4", "N5"],
                       ["C4", "N5", "C6"],
                       ["N5", "C6", "N1"],
                        ["O4'", "C1'", "N1"],
                        ["C3", "C4", "O4"],
                         ["O4", "C4", "N5"],
                         ["C3", "C4", "N4"],
                         ["N4", "C4", "N5"],
                        ["N5", "C6", "O6"],
                        ["O6", "C6", "N1"]
                                            ]

combined_results = [[], [], [], [], [], [], [], []]

# we will hold list of tuple which is (residue, chain, atom1, atom2, distance)


def getDistanceKey(atom1, atom2):
    # this is to return same key for any order of atom1 & atom2, because the sort()
    key = [atom1.get_name(), atom2.get_name()]
    key.sort()
    return ",".join(key)  # convert list to string

def get_atom_by_pattern(data, pattern):
    result = []
    for x in pattern:
        a = [k for k in data if k.get_name() == x ]
        result.extend(a)
    return result

def get_distance(structure_id, filename):

    if filename.suffix == ".cif":
        structure = MMCIFParser().get_structure(structure_id, filename)
    elif filename.suffix == ".pdb":
        structure = PDBParser().get_structure(structure_id, filename)

    dist_by_type = defaultdict(list)

    for model in structure.get_list():
        for chain in model.get_list():
            residues = chain.get_list()  # .get_list gets a list of the subjects children. the children of chain is residues
            for i in range(len(residues)):  # use len to loop through individual residues by index
                if residues[i].get_id()[0] != " ": # exclude bad residue
                    continue
                # following is to tell if its a purine or pyrimidine
                is_pyrimidine = True
                if residues[i].has_id("O4'") and residues[i].has_id("C1'") and residues[i].has_id("N9") \
                        and residues[i].has_id("C4"):
                    is_pyrimidine = False
                # get distance (residue, chain, atom1, atom2, distance)
                for k in range(len(hinge_atoms_base)):
                    atoms = get_atom_by_pattern(residues[i], hinge_atoms_base[k])
                    if len(atoms) == 3:
                        for j in range(len(atoms)-1):
                            dist_by_type[getDistanceKey(atoms[j + 1], atoms[j])].append((atoms[j + 1] - atoms[j], residues[i].get_id(), chain, structure_id))  # a,b = d
                if is_pyrimidine:
                    for k in range(len(hinge_atoms_py)):
                        atoms = get_atom_by_pattern(residues[i], hinge_atoms_py[k])
                        if len(atoms) == 3:
                            for j in range(len(atoms) - 1):
                                dist_by_type[getDistanceKey(atoms[j + 1], atoms[j])].append((atoms[j + 1] - atoms[j], residues[i].get_id(), chain, structure_id))
                else:
                    for k in range(len(hinge_atoms_pu)):
                        atoms = get_atom_by_pattern(residues[i], hinge_atoms_pu[k])
                        if len(atoms) == 3:
                            for j in range(len(atoms) - 1):
                                dist_by_type[getDistanceKey(atoms[j + 1], atoms[j])].append((atoms[j + 1] - atoms[j], residues[i].get_id(), chain, structure_id))
                # do O3' - P distance
                if i > 0 and residues[i-1].has_id("O3'") and residues[i].has_id("P"):
                    if residues[i].get_id()[1]-residues[i-1].get_id()[1] == 1:
                        O3 = residues[i-1]["O3'"]
                        P = residues[i]["P"]
                        dist_by_type["O3',P"].append((O3 - P, residues[i].get_id(), chain, structure_id))
    # return dictionary of atom-pair -> [values], like O4C1 -> [(1.2, residu2, ..), (1.1, residu3, ..)..]
    # always carry more info as a tuple, so later you can make use of it.
    return dist_by_type

def filter_outliers(data1):
    results = []
    distances = [x[0] for x in data1 ]  # distances is ist of all distance values
    stdev = np.std(distances)
    meanv = np.mean(distances)
    for i in range(len(distances)):
        z_score = (distances[i] - meanv) / stdev
        if (z_score < -2) or (
                z_score > 2):  # if the z score is desrible, we convert the tuple into a list and add the z score value
            data = list(data1[i])
            data.append("{:.2f}".format(z_score))
            data.append(stdev)
            data.append(meanv)
            results.append(data)
    return results

def draw_outliers_table(data2, structure_id):
    colTitle = ["distance value", "atom1, atom2", "residue", "chain", "z-score", "standered deviation", "mean"]
    im = plt.imread("C:/David/Books/Torsion_angles.jpg")
    fig, ax = plt.subplots(2, 1, figsize=(12, 20))  # This will make a table of graphs consisting of two rows and 1 column.
    ax[0].imshow(im)  # put image into this chart
    ax[0].axis('off')
    ax[1].axis('off')  # hide axis, just show table
    ax[1].table(cellText=data2,  colLabels=colTitle, loc='center')
    #plt.show()
    plt.savefig(f"C:/Temp/z_Score/{structure_id}_z_score_distance_values.pdf", format="pdf", bbox_inches="tight")


if __name__ == '__main__':
    listOfFiles = [p for p in pathlib.Path("c:/temp/test").iterdir() if p.is_file() and (p.suffix == ".cif" or p.suffix == ".pdb")]
    listOfIds = [f.stem for f in listOfFiles]  # .stem return filename without .extension and path
    for i in range(len(listOfFiles)):
        data1 = get_distance(listOfIds[i], listOfFiles[i])
        #data2 = filter_outliers(data1)
        #draw_outliers_table(data2, listOfIds[i])