#!/usr/bin/env python3
import re
import math
import argparse
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Select best K from NGSadmixture. ")
    parser.add_argument("--inf", required=True, help="Path of the input folder with all the NGSadmixture log files, including the run name, with format /PATH/TO/FODELR/RUNNAME, and the files should have names like RUNNAME_K_R.log, where NAME can be anything K is the K value tested and R is the repetition run for that K.")
    parser.add_argument("--outf", type=str, default="DeltaK", help="Output file name, default DeltaK")
    return parser.parse_args()

args = parse_args()
inputf = args.inf 
outputf = args.outf

NGS_logs = glob.glob(f"{inputf}*.log")
Kmax=0
Rmax=0
for logF in NGS_logs:
    K, R, = re.split(r'/|_|\.', logF )[-3:-1]
    K = int(K)
    R = int (R)
    if K > Kmax:
        Kmax = K
    if R > Rmax:
        Rmax =R

NadPL=[]
NadST=[]
for K in range(1,Kmax+1):
    Ksv=[]
    for R in range(1, Rmax+1):
        inf=(f'{inputf}_{K}_{R}.log')
        with open(inf, 'r') as inf:
            for line in inf:
                if "best" in line:
                    part = re.split(r'[ =]', line)
                    Ksv.append(float(part[2]))
    Kave = np.mean(Ksv)
    Ksd = np.std(Ksv)
    NadPL.append(float(Kave))
    NadST.append(float(Ksd))

dK=[]
if Rmax < 2:
    for i in range(0,Kmax):
        myK = abs((max(NadPL)-NadPL[0])*(i+1-1) - (len(NadPL)-1)*(NadPL[i]-NadPL[0])) / math.sqrt(((max(NadPL)-NadPL[0])**2) + (len(NadPL)-1)**2)
        dK.append(myK)

if Rmax > 1 :
    dK = [0] * (Kmax)
    for i in range (1, Kmax-1):
        DPPK = NadPL[i +1]- 2*NadPL[i] + NadPL[i-1]
        if NadST[i] == 0:
            NadST[i] = NadST[i] + 0.1
        DKN = DPPK / NadST[i]
        dK[i]=DKN

with open(outputf, "w") as of:
    for k, delta_k in enumerate(dK):
        of.write(f"K{k + 1}: {delta_k:.6f}\n")


print("ΔK Calculation Results:")
for k, delta_k in enumerate(dK):
    print(f"K{k}: {delta_k}")

optimal_K = np.argmax(dK)+1
print(f"\nOptimal K: {optimal_K}")
