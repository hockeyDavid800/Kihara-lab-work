import find_deviant_distances as distance
import find_deviant_angles as angle
import matplotlib.pyplot as plt
import numpy as np
import pathlib
from matplotlib.backends.backend_pdf import PdfPages
import pdb_to_image as fetch_image
from collections import defaultdict
import csv
import math

# lets say a Angel-info is ((angel_value,"alpha", residues[i].get_id(),residues[i-1].get_id(),chain,structure_id))

# (angel_value,"alpha", residues[i].get_id(),residues[i-1].get_id(),chain,structure_id) is a tuple containing all angle information. These tuples are the index of all_file_angles
all_file_angles = [[],[],[],[],[],[],[],[]]  # each index is an angles properties. the first index is a list of the residues and chains of all alpha angles
all_file_distances = defaultdict(list) # This defines all_file_distances as a  default dictionary. each index of the dictionary is a list. A default dictionary is a dictionary where
                                        # a defualt dictionary will create a key for any non defined value I add
angle_type =["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "pyrimidine", "purine"]

def import_csv(file):  #this function reads in the standered deviation and mean values for each angle/ distance type in order to comparte them to a pdb's indidiual
                       # values and calculate teh z score
    csv_dict = {} #first i create an empty dictionary
    with open(file, 'r') as csvfile:
        reader = csv.DictReader(csvfile) #reader is all rows in the csv file that can now be operated on by python
        for row in reader:
            csv_dict[row['type']] = (float(row['standard deviation']), float(row['mean']))  # we copy the entire csv into a python dictionary so we can calculate z score earlier
    return csv_dict

def getAngelDiff2(anglea, angleb):
    # 10 - 50, same results as 10 - 330
    anglea = ( anglea * math.pi) / 180
    angleb = ( angleb * math.pi) / 180
    a_x = np.cos(anglea)
    a_y = np.sin(anglea)
    b_x = np.cos(angleb)
    b_y = np.sin(angleb)
    angel_diff = (np.arctan2(a_y-b_y, a_x-b_x) * 180) / math.pi
    return angel_diff

def getAngelDiff(anglea, angleb):
    # 10 - 50, same results as 10 - 330
    anglea = ( anglea * math.pi) / 180
    angleb = ( angleb * math.pi) / 180
    #atan2(sin(x - y), cos(x - y))
    angel_diff = abs(np.arctan2(np.sin(anglea-angleb), np.cos(anglea-angleb)))
    if angel_diff > math.pi:
        angel_diff = 2 * math.pi - angel_diff
    return (angel_diff *180)/ math.pi


def getMeanStd_angle():
    """
    all_file_outliners is list of (angle, type, residue, chain, structure-id, z-score, stdev, meanv)
    :return:
    """
    all_file_outliners = []
    angle_dict = import_csv("C:/Temp/mean_and_stdev/angles.csv")  #angle_dict is the standered deviation and mean for all angle types
    for i in range(len(all_file_angles)): # i is each angle type. ex alpha, beta
        angles = [x[0] for x in all_file_angles[i]]  # x is a tuple of all information regarding each angle type. Therefore x[0] is angle's numerical value.
        (stdev, meanv) = angle_dict[angle_type[i]] # this extracts the stdev and meanv columns for each angle type
        for j in range(len(all_file_angles[i])): #for each angle type
            z_score = getAngelDiff(angles[j], meanv) /stdev
            #z_score = (angles[j] - meanv) / stdev   # get the z score for all alpha angles
            if (z_score < -2) or (z_score > 2):   # if an z score of an angle value is < -2 or >2
                data = list(all_file_angles[i][j])    # all_file_angles[0] is the entire angel-info for all alpha angles
                                                      # therefore, all_file_angles[0][0] is the angle-info for the first alpha angle.
                                                       #overall, all_file_angles[i][j] is the Angel-info for each individual angle
                                                      # we convert all_file_angles into a list because it is a tuple which can't be changed
                data[0] = "{:.2f}".format(data[0])   # we convert all angle values into having two decimal digits
                data.append("{:.2f}".format(z_score))
                data.append("{:.3f}".format(stdev))
                data.append("{:.3f}".format(meanv))
                all_file_outliners.append(data)
    return all_file_outliners

def getMean_distance():
    all_file_outliners = []
    distance_dict = import_csv("C:/Temp/mean_and_stdev/distances.csv")
    for k, v in all_file_distances.items():  # k is the key. it represents the atom pair. v is the list of attributes for that atom pair ex: distance, stdev, etc
        values = [x[0] for x in v]
        (stdev, meanv) = distance_dict[k]
        #print(stdev, meanv)

        for j in range(len(v)):
            z_score = (values[j] - meanv) / stdev  # get the z score for all alpha angles
            if (z_score < -2) or (z_score > 2):  # if an z score of an angle value is < -2 or >2
                data = list(v[j]) # convert tuple into list, so you can modify/append element in list
                data[0] = "{:.2f}".format(data[0])
                data.append(k)
                data.append("{:.2f}".format(z_score))
                data.append("{:.3f}".format(stdev))
                data.append("{:.3f}".format(meanv))
                all_file_outliners.append(data)
    return all_file_outliners



if __name__ == '__main__':
    listOfFiles = [p for p in pathlib.Path("c:/temp/david").iterdir() if p.is_file() and (p.suffix == ".cif" or p.suffix == ".pdb")]  #defines
                                                                                                                                     # listofFiles as all pdb files
    #                                                                                                                               in the selected folder
    listOfIds = [f.stem for f in listOfFiles]  # .stem returns the filename without .pdb
    #pdb_list = zip(listOfFiles,listOfIds)
    #fetch_image.makeImage(pdb_list)

    for i in range(len(listOfFiles)):
        data1 = distance.get_distance(listOfIds[i], listOfFiles[i]) #data1 contains all atom distance values per pdb file
        for k, v in data1.items(): #the key is an atom pair like P1,P2. The value is the atom pair's attributes, which are the distance,id, chain, structure id
            all_file_distances[k].extend(data1[k]) #this updates the all_file_distances list to include the atom pairs and their attributes. it matches the old and new values by the key(atom pair)
        data3 = angle.torsion_angles(listOfIds[i], listOfFiles[i])  # this is all torsion angle values per torison type for each pdb file
        for j in range(len(all_file_angles)):
            all_file_angles[j].extend(data3[j])  #this updates the empty all_file_angles list to include all the torsion types and their angles

    angle_outliers = getMeanStd_angle() # angle_outliers is the list of mean and standered devation angles after filtering them from the list of all angles
    distance_outliers = getMean_distance() #distance_outliers is the list of mean and stdev distance values after filtering them from the list of all distances

    for i in range(len(listOfIds)):
        with PdfPages(f"C:/Temp/PDB_overview/{listOfIds[i]}_overview.pdf") as pdf:
            im = plt.imread(f"{listOfFiles[i]}.png")
            plt.figure(figsize=(6, 8))
            plt.imshow(im)  # put image into this chart
            plt.axis('off')  # hide axis, just show table
            plt.tight_layout()
            plt.title('Page One')
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

            df = [x for x in angle_outliers if x[5] == listOfIds[i]] #makes sure the structure ids match
            if len(df) > 0:
                angleTitle = ["angle value", "angle type", "residue 1", "residue 2", "chain", "structure id",
                              "z -score",
                              "standered deviation", "mean"]
                plt.figure(figsize=(16, 50))
                plt.axis('off')
                table = plt.table(cellText=df, colLabels=angleTitle, loc='center')
                table.scale(1.2, 1.2)
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()


            df = [x for x in distance_outliers if x[3] == listOfIds[i]]
            if len(df) > 0:
                distanceTitle = ["distance value", "residue", "chain", "structure id", "atom1, atom2", "z-score",
                                 "standered deviation", "mean"]
                plt.figure(figsize=(16, 20))
                plt.axis('off')  # hide axis, just show table
                table2 = plt.table(cellText=df, colLabels=distanceTitle, loc='center')
                table2.scale(1.2, 1.2)
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()








