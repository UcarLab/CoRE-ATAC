import numpy as np
import pandas as pd
import argparse
import os

wd = os.getcwd()

parser = argparse.ArgumentParser(description='Formats a peak file for CoRE-ATAC')
parser.add_argument('peakfile', type=str, help='Peak file to be formatted.')
parser.add_argument('fastachromosomes', type=str, help='File listing fasta chromosomes.')
parser.add_argument('outputfile', type=str, help='The output file.')

parser.add_argument('-e', dest='extend', type=int, default=300, help='Amount to extend from the peak center')
parser.add_argument('-a', dest='annotation', type=int, default=None, help='Index column for annotations. (Only used if the peaks will be used for training a model)')

args = parser.parse_args()

peaks = pd.read_csv(args.peakfile, sep="\t", header=None).values
extend = args.extend
aindex = args.annotation
outfile = args.outputfile
fachromosomes = pd.read_csv(args.fastachromosomes, sep="\t", header=None).values[:,0]

rv = []
rv2 = []    
for i in range(len(peaks)):
    curpeak = peaks[i]
    curchr = curpeak[0]
    if ~np.isin(curchr, fachromosomes):
        continue

    start = curpeak[1]
    end = curpeak[2]
    center = start+(end-start)/2
    dlstart = center-extend
    dlend = center+extend-1
    if aindex:
        rv.append([curchr, dlstart, dlend, i, start,end, curpeak[aindex]])
        rv2.append([curchr, start,end, curpeak[aindex]])
    else:
        rv.append([curchr, dlstart, dlend, i, start,end])
        rv2.append([curchr, start,end, i])

rv = np.array(rv)
rv2 = np.array(rv2)


pd.DataFrame(rv).to_csv(outfile, sep="\t", index=None, header=None)
pd.DataFrame(rv2).to_csv(outfile+"_original.txt", sep="\t", index=None, header=None)
