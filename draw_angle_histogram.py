import pathlib
import find_deviant_angles as angle
import numpy as np
import math
import csv
import matplotlib.pyplot as plt


RSCB = [[],[],[],[],[],[],[],[]]  # combined all pdb's angle into each type ex: first index is alpha, second is beta, etc
Lab = [[],[],[],[],[],[],[],[]]
angle_type =["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "pyrimidine", "purine"]


def import_csv(file):
    csv_dict = {}  # first i create an empty dictionary
    with open(file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)  # reader is all rows in the csv file that can now be operated on by python
        for row in reader:
            csv_dict[row['type']] = (float(row['standard deviation']), float(
                row['mean']))  # we copy the entire csv into a python dictionary so we can calculate z score earlier
    return csv_dict



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
    for i in range(len(Lab)): # i is each angle type. ex alpha, beta
        angles = [x[0] for x in Lab[i]]  # x is a tuple of all information regarding each angle type. Therefore x[0] is angle's numerical value.
        (stdev, meanv) = angle_dict[angle_type[i]] # this extracts the stdev and meanv columns for each angle type
        for j in range(len(Lab[i])): #for each angle type
            z_score = getAngelDiff(angles[j], meanv) /stdev
            #z_score = (angles[j] - meanv) / stdev   # get the z score for all alpha angles
            if (z_score < -2) or (z_score > 2):   # if an z score of an angle value is < -2 or >2
                data = list(Lab[i][j])    # all_file_angles[0] is the entire angel-info for all alpha angles
                                                      # therefore, all_file_angles[0][0] is the angle-info for the first alpha angle.
                                                       #overall, all_file_angles[i][j] is the Angel-info for each individual angle
                                                      # we convert all_file_angles into a list because it is a tuple which can't be changed
                #data[0] = data[0]   # we convert all angle values into having two decimal digits
                data.append("{:.2f}".format(z_score))
                data.append("{:.3f}".format(stdev))
                data.append("{:.3f}".format(meanv))
                all_file_outliners.append(data)
    return all_file_outliners



if __name__ == '__main__':
    listOfFiles = [p for p in pathlib.Path("c:/temp/RSCB").iterdir() if p.is_file() and (p.suffix == ".cif" or p.suffix == ".pdb")]
    listOfIds = [f.stem for f in listOfFiles]  # .stem return filename without .extension and path
    for i in range(len(listOfFiles)):  # for each file
        data3 = angle.torsion_angles(listOfIds[i], listOfFiles[i])  # this is all torsion angle values
        for j in range(len(RSCB)):
            RSCB[j].extend(data3[j])

    listOfLab = [p for p in pathlib.Path("c:/temp/david").iterdir() if p.is_file() and (p.suffix == ".cif" or p.suffix == ".pdb")]
    listOflab_ids = [f.stem for f in listOfLab]  # .stem returns the filename without .pdb
    for i in range(len(listOfLab)):
        data2 = angle.torsion_angles(listOflab_ids, listOfLab[i])  # this is all torsion angle values per torison type for each pdb file
        for j in range(len(Lab)):
            Lab[j].extend(data2[j])  #this updates the empty all_file_angles list to include all the torsion types and their angles

    angle_outliers = getMeanStd_angle()

    plt.figure(figsize=(16, 12))
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.tight_layout()  # make sure the titles don't overlap
    plt.subplots_adjust(hspace=0.4, wspace=0.3)  # space between subplots

    for i in range(len(Lab)):
        plt.subplot(4, 2, i + 1)
        base = [x[0] for x in RSCB[i]]
        outliers = [x[0] for x in angle_outliers if x[1] == angle_type[i] ]
        normal = [x[0] for x in Lab[i] if x[0] not in outliers]

        bins = np.linspace(-180, 180, 100)
        plt.hist(base, bins, alpha=0.5, label='RSCB angles')
        plt.hist(outliers, bins, alpha=0.5, label='outlier angles')
        plt.hist(normal, bins, alpha=0.5, label='lab data')
        plt.legend(loc='upper right')
        plt.title(f'comparison distribution of {angle_type[i]} ', fontsize=10)

    plt.savefig(f"C:/Temp/check_work/compare_distribution2.pdf", format="pdf", bbox_inches="tight")



