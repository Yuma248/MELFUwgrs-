#!/usr/bin/python

#This script is from https://github.com/jarobin/CaliforniaCondorGenome2021/blob/main/SlidingWindowHet.py with some samll modifications 


import sys
import pysam
import os
import gzip
import argparse

parser = argparse.ArgumentParser(description='This script counts number of called genotypes and number of heterozygotes per sample in sliding windows',epilog=("Example: python ./SWHet.py -v input.vcf.gz -c Chromolenght -w 100000 -s 10000"))
parser.add_argument("-v","--IVCF", type=str,help="Input VCF, it can be single- or multi-sample but should be filtered and it should have extension .vcf or .vcf.gz")
parser.add_argument("-c", "--chrlen", type=str, help="Chromosome lengths information, tab delimited file with two columns, first the name of chromosomes as they are in the VCF, and then the size of each of them")
parser.add_argument("-w", "--widsiz", type=int, default=1000000, help="window size, to calculate Heterozygosity, default 1000000")
parser.add_argument("-s", "--stesiz", type=int, default=1000000, help="step size, the number of bp that the window move to calculate heterozygosity, if this is equal to the window size, the windows will not overlap, default 1000000")

args = parser.parse_args()

if args.IVCF and args.chrlen:
# Open input file and make sure the VCF file is indexed (if not, create index)
	filename = args.IVCF

	if  os.path.exists(filename):
		if filename.endswith("gz") :
			VCF = gzip.open(filename, 'r')

		elif filename.endswith("vcf") : 
			VCF = open(filename, 'r')

		else :
			parser.error("VCF file is not in the righ format. ")
	else:
		parser.error("Input  VCF does not exist. Use -h argument to check each argument. ")


# Generate a dictionary with chromosomes and chromosome lengths
	chrfile = args.chrlen
	if  os.path.exists(chrfile):
		with open(chrfile, 'r') as cc:
			chrom_size={line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in cc}
			CHRL = []
			CHRL =chrom_size.keys()
		cc.close()
	else:
		parser.error("Input  chromosome file does not exist. Use -h argument to check each argument. ")


	if not os.path.exists("%s.tbi" % filename):
		pysam.tabix_index(filename, preset="vcf")
		parsevcf = pysam.Tabixfile(filename)
	else:
		parsevcf = pysam.Tabixfile(filename)


# Set variables
	window_size = int(args.widsiz)
	step_size = int(args.stesiz)





# Get list of samples from VCF file header
	samples=[]
	for line in VCF:
		mline = line.decode('UTF-8')
		if mline.startswith('##'):
			pass
		else:
			for i in mline.split()[9:]:
				samples.append(i)
			break

# Create output file
	outname = filename.replace(".vcf","").replace(".gz","").replace(".recode","")
	output = open(outname + '_het_%swin_%sstep.txt' % (window_size, step_size), 'w')
	output.write('chrom\twindow_start\tsites_total\tcalls_%s\thets_%s\thetp_%s\n' % ('\tcalls_'.join(samples), '\thets_'.join(samples), '\thetp_'.join(samples)))

# Fetch a region, ignore sites that fail filters, tally genotype calls and heterozygotes        
	def snp_cal(chrom,window_start,window_end):
		try:
			rows = tuple(parsevcf.fetch(region="%s:%s-%s" % (chrom, window_start, window_end), parser=pysam.asTuple()))    
		except  Exception as e:
			return
		sites_total=0
		calls=[0]*len(samples)
		hets=[0]*len(samples)
		hetsamp=[0]*len(samples)
		for line in rows:
			sites_total+=1
			for i in range(0,len(samples)):
				if line[i+9][:1]=='.': continue
				calls[i]+=1
				GT=line[i+9].split(':')[0]
				if '/' in GT: sp='/'
				if '|' in GT: sp='|'
				if GT.split(sp)[0]!=GT.split(sp)[1]: hets[i]+=1
		for i in range(0,len(samples)):
			if (hets[i] > 0): hetsamp[i]=hets[i]/calls[1]
			else: hetsamp[i]=0

		output.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom,window_start,sites_total,'\t'.join(map(str,calls)),'\t'.join(map(str,hets)),'\t'.join(map(str,hetsamp)) ) )



	for chromname in CHRL:
		chrom = chromname.rstrip('\n')

# Get start and end positions of chromosome
		for mline in VCF:
			line = mline.decode('UTF-8')
			if line[0] != '#':
#				start_pos = int(line.strip().split()[1])
				start_pos = 1
				end_pos = int(chrom_size[chrom])
				#print (chrom,start_pos, end_pos)
				#print (end_pos)
				break
# Initialize window start and end coordinates
		window_start = start_pos
		window_end = start_pos+window_size-1


# Calculate stats for window, update window start and end positions, 
# repeat to end of chromosome
		while window_end <= end_pos:
			if window_end < end_pos:
				snp_cal(chrom,window_start,window_end)
				window_start = window_start + step_size
				window_end = window_start + window_size - 1
			else:
				snp_cal(chrom,window_start,window_end)
				break
		else:
			window_end = end_pos
			snp_cal(chrom,window_start,window_end)







# Close files and exit
	VCF.close()
	output.close()

	exit()


else:
	parser.error("Both, a VCF and a chromosome length files are required, run the script with the -h argument to check each argument. ")




