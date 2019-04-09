import os
from os import listdir
from os.path import isfile, join

import sys, time
import argparse

from dataset import * #helper class to hanlde the dataset

import pandas as pd


def npDoubletsLoad(path):
    print ("======================================================================")

    start = time.time()
    bal_dir = path + "/bal_data/"
    new_dir = path + "/original/"

    datafiles = [f for f in listdir(path) if (isfile(join(path, f)) and  f.lower().endswith(("txt","gz")) and "dnn_doublets" in f)]

    print("Loading " + str(len(datafiles)) + " dataset file(s) . . .")

    print("Saving  original in   : " + new_dir)
    print("Saving  balanced in     : " + bal_dir)

    if not os.path.exists(bal_dir):
        os.makedirs(bal_dir)
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    idName = ""

    for p in path.split("/"):
        if "runs" in p:
            idName = p


    print(idName)

    listdata = []
    for no,d in enumerate(datafiles):
        if os.stat(path + "/" + d).st_size == 0:
                print("File no." + str(no+1) + " " + d + " empty.Skipping.")
                continue
        with open(path + "/" + d, 'rb') as df:
            print("Reading file no." + str(no+1) + ": " + d)
            dfDoublets = pd.read_table(df, sep="\t", header = None)

            print("--Dumping unbalanced data")
            dfDoublets.columns = dataLab
            dfDoublets.to_hdf(new_dir + idName + "_" + d.replace(".txt",".h5"),'data',append=False,complib="bzip2",complevel=9)

            ##balanceddata
            print("--Dumping balanced data")
            theData = Dataset([])
            theData.from_dataframe(dfDoublets)
            theData.balance_data()
            theData.save(bal_dir + idName + "_bal_" + d.replace(".txt",".h5"))

            os.remove(path + "/" + d)


    end = time.time()
    print ("======================================================================")
    print ("\n - Timing : " + str(end-start))



if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog="dataToHDF")
    parser.add_argument('--read', type=str, default="doublets/",help='files path')
    args = parser.parse_args()

    npDoubletsLoad(args.read)
