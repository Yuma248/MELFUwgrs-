#!/usr/bin/python

import sys
import os
import gzip
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='This script create take a subset of individual from a beagle file.',epilog=("Example: python IndBeagle.py -b input.beagle.gz -l bamlist -i indlist -o output"))
parser.add_argument("-b","--beagle", type=str,help="Input beagle file, a multi-sample beagle files it could be  .gz compressed")
parser.add_argument("-l", "--bamlist", type=str, help="The bamlist file used to produce the beagle file, tab delimited file with two columns, first the name of chromosomes as they are in the VCF, and then the size of each of them")
parser.add_argument("-i", "--indlist", type=str, help="A list of samples that you want to subsample from the beagle file, the names should be in the same format as in the bamfile excluding paths or folders")
parser.add_argument("-o", "--output", type=str, default="NewBeagle.gz", help="Path and name to the new beagle file with a .gz extension")

args = parser.parse_args()

if args.beagle and args.bamlist and args.indlist:
# Open input file
	filename = args.beagle
	if  os.path.exists(filename):
		if filename.endswith("gz") :
			with gzip.open(filename, 'r') as file:
				BEA = pd.read_csv(file, sep='\t')

		elif filename.endswith("beagle") : 
			with open(filename, 'r') as file:
				BEA = pd.read_csv(file, sep='\t')

		else :
			parser.error("Can no open beagle file. ")

	else:
		parser.error("Beagle file does not exist. Use -h to check argument. ")


# Generate a dictionary of individual 
	indind = args.bamlist
	if  os.path.exists(indind):
		with open(indind, 'r') as idi:
			lni= 0
			ind_ind={}
			for line in idi:
				lni += 1
				ind_ind[line.strip().split('/')[-1]] = lni
			IND = []
			IND =ind_ind.keys()
		idi.close()
	else:
		parser.error("Input bamlist does not exist. Use -h to check each argument. ")

# Get list of samples to subsaple form indlist
	samples=[]
	inds = args.indlist
	if  os.path.exists(inds):
		with open(inds, 'r') as idl:
			for line in idl:
				samples.append(line.strip())
		idl.close()
	else:
		parser.error("Input indlist does not exist. Use -h to check each argument. ")





#Select the index to be extracted
	mycol=[]
	mycol.extend(range(0,3))
	for i in samples:
		mycol.extend(range(ind_ind[i]*3,ind_ind[i]*3+3))

	newB = BEA.iloc[:, mycol]

# Create output file
	outname = args.output
	newB.to_csv(outname, sep ='\t',  index=False, compression='gzip')




	exit()


else:
	parser.error("A beagle, bamlist and a subsample list files are required, run the script with the -h to check each argument.")


