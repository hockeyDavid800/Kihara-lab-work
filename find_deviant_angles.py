import numpy
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.vectors import Vector, calc_dihedral, calc_angle
import matplotlib.pyplot as plt
import numpy as np
import pathlib

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

combined_results = [[], [], [], [], [], [], [], []]

def getTorsion(aList): #calculates the torsion angles for each set of 4 angles
    vectors = [x.get_vector() for x in aList]
    angle = calc_dihedral(*vectors)
    return angle * 180 / 3.1416

def torsion_angles(structure_id, filename):

    if filename.suffix == ".cif":
        structure = MMCIFParser().get_structure(structure_id, filename)
    elif filename.suffix == ".pdb":
        structure = PDBParser().get_structure(structure_id, filename)

    angle_value_list = [[], [], [], [], [], [], [], []]

    for model in structure.get_list():
        for chain in model.get_list():
            residues = chain.get_list()  # .get_list gets a list of the subjects children. the children of chain is residues
            for i in range(len(residues)):  # i is each residue
                if residues[i].get_id()[0] != " ": # exclude bad residue
                    continue
                # get the alpha angle
                if i > 0 and residues[i - 1].has_id("O3'") and residues[i].has_id("P") and residues[i].has_id("O5'") \
                        and residues[i].has_id("C5'") and residues[i - 1].get_id()[1] + 1 == residues[i].get_id()[1]:
                    if residues[i].get_id()[1] - residues[i - 1].get_id()[1] == 1:
                        d = [residues[i - 1]["O3'"], residues[i]["P"], residues[i]["O5'"], residues[i]["C5'"]]
                        a = getTorsion(d)
                        angle_value_list[0].append((a,"alpha", residues[i].get_id(),residues[i-1].get_id(),chain,structure_id))
                # get beta angle vale, all in same residue, no index issue, no seq issue
                if residues[i].has_id("P") and residues[i].has_id("O5'") and residues[i].has_id("C5'") \
                        and residues[i].has_id("C4'"):
                    d = [residues[i]["P"], residues[i]["O5'"], residues[i]["C5'"], residues[i]["C4'"]]
                    a = getTorsion(d)
                    angle_value_list[1].append((a,"beta", residues[i].get_id(), residues[i].get_id(),chain, structure_id))
                # get gamma
                if residues[i].has_id("O5'") and residues[i].has_id("C5'") and residues[i].has_id("C4'") \
                        and residues[i].has_id("C3'"):
                    d = [residues[i]["O5'"], residues[i]["C5'"], residues[i]["C4'"], residues[i]["C3'"]]
                    a = getTorsion(d)
                    angle_value_list[2].append((a,"gamma", residues[i].get_id(),residues[i].get_id(),chain, structure_id))
                # get delta
                if residues[i].has_id("C5'") and residues[i].has_id("C4'") and residues[i].has_id("C3'") \
                        and residues[i].has_id("O3'"):
                    d = [residues[i]["C5'"], residues[i]["C4'"], residues[i]["C3'"], residues[i]["O3'"]]
                    a = getTorsion(d)
                    angle_value_list[3].append((a,"delta", residues[i].get_id(),residues[i].get_id(),chain, structure_id))
                # get epsilon
                if i < len(residues) - 1 and residues[i].has_id("C4'") and residues[i].has_id("C3'") \
                        and residues[i].has_id("O3'") and residues[i + 1].has_id("P") \
                        and residues[i].get_id()[1] + 1 == residues[i + 1].get_id()[1]:
                    if residues[i+1].get_id()[1] - residues[i].get_id()[1] == 1:
                        d = [residues[i]["C4'"], residues[i]["C3'"], residues[i]["O3'"], residues[i + 1]["P"]]
                        a = getTorsion(d)
                        angle_value_list[4].append((a,"epsilon", residues[i].get_id(),residues[i+1].get_id(),chain, structure_id))
                # get zeta
                if i < len(residues) - 1 and residues[i].has_id("C3'") and residues[i].has_id("O3'") \
                        and residues[i + 1].has_id("P") and residues[i + 1].has_id("O5'") \
                        and residues[i].get_id()[1] + 1 == residues[i + 1].get_id()[1]:
                    if residues[i+1].get_id()[1] - residues[i].get_id()[1] == 1:
                        d = [residues[i]["C3'"], residues[i]["O3'"], residues[i + 1]["P"], residues[i + 1]["O5'"]]
                        a = getTorsion(d)
                        angle_value_list[5].append((a,"zeta", residues[i].get_id(),residues[i+1].get_id(),chain, structure_id))

                is_pyrimidine = True
                #get chi purine
                if residues[i].has_id("O4'") and residues[i].has_id("C1'") and residues[i].has_id("N9") \
                        and residues[i].has_id("C4"):
                    is_pyrimidine = False
                    d = [residues[i]["O4'"], residues[i]["C1'"], residues[i]["N9"], residues[i]["C4"]]
                    a = getTorsion(d)
                    angle_value_list[7].append((a,"chi purine", residues[i].get_id(),residues[i].get_id(),chain, structure_id))
                # get chi pyrimidine
                if is_pyrimidine and residues[i].has_id("O4'") and residues[i].has_id("C1'") \
                        and residues[i].has_id("N1") and residues[i].has_id("C2"):
                    d = [residues[i]["O4'"], residues[i]["C1'"], residues[i]["N1"], residues[i]["C2"]]
                    a = getTorsion(d)
                    angle_value_list[6].append((a,"chi pyrimidine", residues[i].get_id(),residues[i].get_id(),chain, structure_id))
    # returned data is list of angel_type, each angle_type is a list of (angle_vale, angle_name, residue, chain, structure-id)
    return angle_value_list


def find_outliers_table(values):
    results = []
    for i in range(len(values)):
        angles = [x[0] for x in values[i]]  # for each tuple, get a value out of it
        stdev = np.std(angles)
        meanv = np.mean(angles)
        # f.write(f" the standered deviation is {stdev}\n")
        # f.write(f" the mean is {meanv}\n")
        for j in range(len(angles)):
            z_score = (angles[j] - meanv) / stdev
            if (z_score < -2) or (z_score > 2):
                data = list(values[i][j])
                data[0] = "{:.2f}".format(data[0])
                data.append("{:.2f}".format(z_score))
                data.append(stdev)
                data.append(meanv)
                results.append(data)
    return results


def draw_outliers_table(data2, structure_id):
    colTitle = ["angle value", "angle type", "residue 1", "residue 2", "chain", "structure id", "z -score", "standered deviation", "mean"]
    im = plt.imread("C:/David/Books/Torsion_angles.jpg")
    fig, ax = plt.subplots(2, 1,
                           figsize=(12, 12))  # This will make a table of graphs consisting of two rows and 1 column.
    ax[0].imshow(im)  # put image into this chart
    ax[0].axis('off')
    ax[1].axis('off')  # hide axis, just show table
    ax[1].table(cellText=data2, colLabels=colTitle, loc='center')
    plt.tight_layout()
    # plt.show()
    plt.savefig(f"C:/Temp/z_Score/{structure_id}_z_score_values.pdf", format="pdf", bbox_inches="tight")


if __name__ == '__main__':
    listOfFiles = [p for p in pathlib.Path("c:/temp/test").iterdir() if p.is_file() and (p.suffix == ".cif" or p.suffix == ".pdb")]
    listOfIds = [f.stem for f in listOfFiles]  # .stem return filename without .extension and path
    for i in range(len(listOfFiles)): #for each pdb file
        values = torsion_angles(listOfIds[i], listOfFiles[i]) #finds the torsion angles for each pdb file
        data2 = find_outliers_table(values)
        draw_outliers_table(data2, listOfIds[i])