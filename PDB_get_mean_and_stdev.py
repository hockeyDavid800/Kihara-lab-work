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

all_file_angles = [[],[],[],[],[],[],[],[]]  # combined all pdb's angle into each type ex: first index is alpha, second is beta, etc
angle_type =["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "pyrimidine", "purine"]
all_file_distances = defaultdict(list)


def getMeanAngel(angelList): # this function takes the angel values and calculates their mean for each torsion type by using calculus
    x_list = []
    y_list = []
    for a in angelList:
        r = (a * math.pi) / 180  # converted to radius
        x_list.append(np.cos(r))
        y_list.append(np.sin(r))
    x_mean = np.mean(x_list)
    y_mean = np.mean(y_list)
    mean_angle = np.arctan2(y_mean, x_mean) * 180 / math.pi
    stdev = math.sqrt(-1 * math.log(y_mean * y_mean + x_mean * x_mean)) * 180 / math.pi
    ##  standard deviation of angles = sqrt(-log R²)
    return ("{:.2f}".format(mean_angle), "{:.2f}".format(stdev))

def makeMeanStd_angle():
    """
    all_file_outliners is list of (angle, type, residue, chain, structure-id, z-score, stdev, meanv)
    :return:
    """
    with open('c:/temp/mean_and_stdev/angles.csv', 'w', newline='', encoding='UTF8') as f:
        writer = csv.writer(f)
        header = ["type", "standard deviation", "mean"]
        writer.writerow(header)
        for i in range(len(all_file_angles)):  # i is each angle type. ex alpha, beta
            angles = [x[0] for x in all_file_angles[i]]  # for each tuple, get a value out of it
            (meanv, stdev) = getMeanAngel(angles)
            writer.writerow([angle_type[i], stdev, meanv])


def makeMean_distance():
    with open('c:/temp/mean_and_stdev/distances.csv', 'w', newline='', encoding='UTF8') as f:
        writer = csv.writer(f)
        header = ["type", "standard deviation", "mean"]
        writer.writerow(header)
        for k, v in all_file_distances.items():  # k is the key. it represents the atom pair. v is the list of attributes for that atom pair ex: distance, stdev, etc
            values = [x[0] for x in v]

            stdev = np.std(values)
            meanv = np.mean(values)
            writer.writerow([k, stdev, meanv])


if __name__ == '__main__':
    listOfFiles = [p for p in pathlib.Path("c:/temp/RSCB").iterdir() if p.is_file() and (p.suffix == ".cif" or p.suffix == ".pdb")]
    listOfIds = [f.stem for f in listOfFiles]  # .stem return filename without .extension and path

    for i in range(len(listOfFiles)):  # for each file
        data1 = distance.get_distance(listOfIds[i], listOfFiles[i]) #this uses the get_distance function from find_deviant_distances.py
                                                                    # to get the distance values for each atom pair for each inviddual pdb file
        for k, v in data1.items():  #The key is the atom pair. The value is the atom distance value, the chain, residue, and file id
            all_file_distances[k].extend(data1[k]) #updates all_file_distances to contain the characteristics for each atom pair
        data3 = angle.torsion_angles(listOfIds[i], listOfFiles[i])  # this is all torsion angle values
        for j in range(len(all_file_angles)):
            all_file_angles[j].extend(data3[j])

    makeMeanStd_angle()
    makeMean_distance()








