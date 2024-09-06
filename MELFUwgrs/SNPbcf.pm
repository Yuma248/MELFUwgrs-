package SNPbcf;
my $LEVEL = 1;
sub callSNPs2 {

use File::Path qw( make_path );
my @arg = @_;
foreach $ar (@arg){
        if ($ar=~ /^-i/){our $inputfolder= (split(/ /,$ar))[1];}
        elsif ($ar=~ /^-o/){our $outputfolder= (split(/ /,$ar))[1];}
        elsif ($ar=~ /^-nc/){our $nc=(split(/ /,$ar))[1];}
        elsif ($ar=~ /^-rg/){$refgenome=(split(/ /,$ar))[1];}
        elsif ($ar=~ /^-chrf/){$chrf=(split(/ /,$ar))[1];}
}
if (not defined ($inputfolder && $outputfolder && $refgenome)){print "\nThis script will create the bcftools command to call and genotype SNPs for several samples in parallel. It needs the mapped bam files in a folder,and an indexed reference genome. It is recommended to use  just chromosomes or the biggest scaffolds, so you can use a file with the names of chromosomes or scaffolds to be used, default will use the any contig or scaffold with chr in the name.\n\nUsage:\nMELFUwgrs -stp snpcalling2\n\t-i <input folder with mapped bam files with realigned indels>\n\t-o <output folder to save vcf files>\n\t-rg <reference genome>\n\t-snc <number the cores to be used in parallel, recommend to use the number of Chromosomes, default 24>\n\t-chrf <a file with the list of scaffolds to use, scaffold or chromosomes names as they appear in the reference, one per row. Default N and will use the first 30 scaffold in the reference fasta file>\n\nExample:\nMELFUwgrs.pl -stp snpcalling -i ./Yuma/indelrealigned/ -o ./Yuma/rawsnp/ -rg ./Yuma/genome/reference_genome.fasta -snc 23 -chrf chromosomefile\n\n"; exit;}
$snc //= 20;
$chrf //= "N";

foreach my $file (@files){ next unless ($file =~ /\.bam$/); my @fileinf = split (/\./, $file);  push @samplesnames, $fileinf[0];}
our $bamlist = "$inputfolder/bamlist";
if ( !-e $bamlist){
foreach $samplename (@samplesnames){
`echo "$inputfolder\/$samplename\.bam " >> $bamlist`;
}
}
if ($chrf eq "N"){
our @scaffr=`grep \">\" $refgenome | head -n 30 | perl -p -e \'s/>//\' | perl -p -e \'s/\\s\.*\\n/\\n/\'`;
} else {
our @scaffr=`cat $chrf`;
}
our @cores = (1 .. (scalar @scaffr));

our $cmd="parallel --link -j $snc bcftools mpileup -Ou -r {1} -f $refgenome -b $bamlist -d 1000 -a \"FORMAT/DP,FORMAT/AD\" \'|\' bcftools call -vmO z -o $outputfolder/call.{2}.vcf.gz ::: @scaffr ::: @cores";
system($cmd)
