#!//usr/bin/env python3

# ---------------------------------------------------------
# Yuma script to obtain a VCF file from an ANGSD beagle file
#
# options:
#  -h, --help            show this help message and exit
#  -b BEAGLE, --beagle BEAGLE
#                        ANGSD beagle file with header
#  -m MAFS, --mafs MAFS  ANGSD maf file
#  -s SAMPLES, --samples SAMPLES
#                        Sample names file
#  -o OUT, --out OUT     Output VCF
#  -p MINPROB, --minprob MINPROB
#                        Min probability, default 0.9
# --------------------------------------------------------
import gzip
import argparse
import math

def open_file(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def genotype_call(probs, min_prob=0.9):
    best = max(range(3), key=lambda x: probs[x])
    best_prob = probs[best]
    if best_prob < min_prob:
        return "./.", 0
    if best == 0:
        gt = "0/0"
    elif best == 1:
        gt = "0/1"
    else:
        gt = "1/1"
    if best_prob >= 1:
        gq = 99
    else:
        gq = min(99, int(-10 * math.log10(1 - best_prob)))
    return gt, gq

parser = argparse.ArgumentParser(
    description="Convert ANGSD Beagle + MAFS to VCF"
)
parser.add_argument("-b", "--beagle", required=True, help="ANGSD beagle file with header")
parser.add_argument("-m", "--mafs", required=True, help="ANGSD maf file")
parser.add_argument("-s", "--samples", required=True, help="Sample names file")
parser.add_argument("-o", "--out", required=True, help="Output VCF")
parser.add_argument("-p", "--minprob",type=float, default=0.9, help="Min probability, defaul 0.9")
args = parser.parse_args()



# -----------------------------
# Read sample names
# -----------------------------
with open(args.samples) as f:
    samples = [
        x.strip()
        for x in f
        if x.strip()
    ]
# -----------------------------
# Read Beagle header
# -----------------------------
beagle = open_file(args.beagle)
beagle_header = next(beagle).strip().split()
# Number of individuals from Beagle
n_individuals = (len(beagle_header) - 3) // 3
if n_individuals != len(samples):
    raise ValueError(
        f"Mismatch: Beagle contains {n_individuals} individuals "
        f"but samples file contains {len(samples)} names"
    )
print(f"Samples loaded: {len(samples)}")
# -----------------------------
# Read MAFS
# -----------------------------
mafs = open_file(args.mafs)
maf_header = next(mafs).strip().split()
chr_col = maf_header.index("chromo")
pos_col = maf_header.index("position")
ref_col = maf_header.index("major")
alt_col = maf_header.index("minor")
# -----------------------------
# Write VCF header
# -----------------------------
vcf = open(args.out, "w")
vcf.write("##fileformat=VCFv4.2\n")
vcf.write(
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
)
vcf.write(
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
)
vcf.write(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
)
vcf.write("\t".join(samples))
vcf.write("\n")
# -----------------------------
# Convert variants
# -----------------------------
for maf_line, beagle_line in zip(mafs, beagle):
    maf = maf_line.strip().split()
    beag = beagle_line.strip().split()
    chrom = maf[chr_col]
    pos = maf[pos_col]
    ref = maf[ref_col]
    alt = maf[alt_col]
    probs = beag[3:]
    genotypes = []
    for i in range(len(samples)):
        p = list(
            map(
                float,
                probs[i*3:(i+1)*3]
            )
        )
        gt, gq = genotype_call(
            p,
            args.minprob
        )
        genotypes.append(
            f"{gt}:{gq}"
        )
    vcf.write(
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:GQ\t"
        + "\t".join(genotypes)
        + "\n"
    )
vcf.close()

print("Done:", args.out)
