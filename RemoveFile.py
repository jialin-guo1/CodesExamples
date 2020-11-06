import argparse
parser = argparse.ArgumentParser(description="reduce files from a list")
parser.add_argument("-i","--inputfiles", dest="inputfiles", default="datasetMC.txt",help="input files")
parser.add_argument("-o","--outputfiles",dest="outputfiles",default="datasetMC_dir.txt",help="output files")
args = parser.parse_args()

import os

print("get file " + str(args.inputfiles))

infile = open(args.inputfiles)
outfile = open(args.outputfiles,"w")
for line in infile:
    for i in line:
        if(i==".") break
    output.write(line)

outfile.Close()
