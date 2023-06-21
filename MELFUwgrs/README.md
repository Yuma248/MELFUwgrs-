# SNPcallPipe
Pipeline to call snps using a reference genome

# Download
        git clone https://github.com/Yuma248/SNPcallingPipe
  
# Dependecies 

Perl Parallel:::Loops

GNU Parallel

Stacks

AdapterRemoval

Bowtie2

BWA

samtools

bcftools

vcftools


The easiest way to install the dependencies is using conda 

conda create --name SNPcallPipe -c conda-forge -c bioconda perl-parallel-loops parallel stacks adapterremoval bowtie2 bwa samtools bcftools vcftools


# Basic usage:

Usage:

SNPcallPipe

        -stp <You need at least determine what steps you want to run>
        
                indref: <Indexs the reference genome with samtools, picard and bowtie2>
                
                demul: <It will use stacks process_rad script, to demultiples samples base on a barcode file>
                
                trim: <It will use AdapterRemoval to trim and filter reads>
                
                aligment: <It will use bowtie2 or bwq to align reads to a referecne genome>
                
                dedup: <This step will sort sam files, cnvert to bam and mask duplicates>
                
                indelrea: <This step will locally realign indels, although this is not recomended any more>
                
                calling: <This step will use bcftool and mpileup to call variant sites SNP/indel>
                
                filtering: <This step will use vcftools to filter SNPs, I recomend to use this automatically to have an idea of youdata, but play wiht the parameters if you have the time>
                


