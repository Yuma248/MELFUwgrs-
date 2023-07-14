package SNPangsd;
my $LEVEL = 1;
sub callSNPs {

use File::Path qw( make_path );
my @arg = @_;
foreach $ar (@arg){
        if ($ar=~ /^-i/){our $inputfolder= (split(/ /,$ar))[1];}
        elsif ($ar=~ /^-o/){our $outputfolder= (split(/ /,$ar))[1];}
        elsif ($ar=~ /^-nc/){our $nc=(split(/ /,$ar))[1];}
        elsif ($ar=~ /^-rg/){$refgenome=(split(/ /,$ar))[1];}
        elsif ($ar=~ /^-chrf/){$chrf=(split(/ /,$ar))[1];}
}

if (not defined ($inputfolder && $outputfolder && $refgenome)){print "\nThis script will create the ANGSD command to call and genotype SNPs for several samples in parallel, it use a likelihood approach recommended for low coverage genome sequencing. Its need the mapped bam files in a folder,and a indexed reference genome. It is recommended to use  just chromosomes or the biggest scaffolds, so you can use a file with the names of chromosomes or scaffolds to be used, default will use the any contig or scaffold with chr in the name.\n\nUsage:\nMELFUwgrs -stp snpcalling\n\t-i <input folder with mapped bam files with realigned indels>\n\t-o <output folder to save vcf files>\n\t-rg <reference genome>\n\t-snc <number the cores to be used in parallel, recommend to use the number of Chromosomes, default 20>\n\nExample:\nMELFUwgrs.pl -stp snpcalling -i ./Yuma/indelrealigned/ -o ./Yuma/rawsnp/ -rg ./Yuma/genome/reference_genome.fasta -snc 23 -chrf\n\n"; exit;}
if (not defined ($snc)){$snc = 20;}
if (not defined ($chrf)){$chrf = "N";}

opendir(DIR, "$inputfolder") or die "Can not open folder, $!\n" ;
our @files = readdir(DIR);
closedir (DIR);
if ( !-d $outputfolder ) {
    make_path $outputfolder or die "Failed to create path: $outputfolder";
}
our @samplesnames=();
our @pops=();

foreach my $file (@files){ next unless ($file =~ /\.bam$/); my @fileinf = split (/\./, $file);  push @samplesnames, $fileinf[0];}
our $bamlist = "$inputfolder/bamlist";
if ( !-e $bamlist){
foreach $samplename (@samplesnames){
`echo "$inputfolder\/$samplename\.bam " >> $bamlist`;
}
}
if ($chrf eq "N"){
our @scaffr=`grep \">\" $refgenome | head -n 20 | perl -p -e \'s/>//\' | perl -p -e \'s/\\s\.*\\n/\\n/\'`;
} else {
our @scaffr=`cat $chrf`;
}
chomp @scaffr;

$cmd1="parallel -j $snc angsd -r {1} -b $bamlist -ref $refgenome -out $outputfolder\/{1} -GL 2 -doMajorMinor 1 -minMapQ 20 -minQ 20 -doMaf 1 -doBCF 1 -SNP_pval 1e-6 -doCounts 1 -minMaf 0.03 -doGeno 3 -doGlf 2 -uniqueonly 1 -remove_bads 1 -only_proper_pairs 0 -C 50 -baq 1 -doPost 1 -postCutoff 0.9 -geno_minDepth 3 -geno_maxDepth 1000 -nThreads 3 -minInd 17 ::: @scaffr";
print "$cmd1\n";
#`echo $cmd1 \>\> ./TestYuma`; 
system ($cmd1);
#print "$count1\t$count2\n";

}

1;







