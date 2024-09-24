package SNPbcf;
my $LEVEL = 1;
sub callSNPs2 {

use File::Path qw( make_path );
#while (@ARGV){
#        $_=shift @ARGV;                                 
#        if ($_=~ /^-i$/){$inputfolder=shift @ARGV;}     
#        elsif ($_=~ /^-o$/){$outputfolder=shift @ARGV;} 
#        elsif ($_=~ /^-rg$/){$refgenome=shift @ARGV;}
#        elsif ($_=~ /^-chrf$/){$chrf=shift @ARGV;}
#        elsif ($_=~ /^-nc$/){$snc=shift @ARGV;}
#}
my @arg = @_;
foreach $ar (@arg){
        if ($ar=~ /^-i/){our $inputfolder= (split(/ /,$ar))[1];}
        elsif ($ar=~ /^-o/){our $outputfolder= (split(/ /,$ar))[1];}
        elsif ($ar=~ /^-nc/){our $snc=(split(/ /,$ar))[1];}
        elsif ($ar=~ /^-rg/){$refgenome=(split(/ /,$ar))[1];}
        elsif ($ar=~ /^-chrf/){$chrf=(split(/ /,$ar))[1];}
}
if (not defined ($inputfolder && $outputfolder && $refgenome)){print "\nThis script will create the bcftools command to call and genotype SNPs for several samples in parallel. It needs the mapped bam files in a folder,and an indexed reference genome. It is recommended to use  just chromosomes or the biggest scaffolds, so you can use a file with the names of chromosomes or scaffolds to be used, default will use the first 24 contig or scaffold in the reference.\n\nUsage:\nMELFUwgrs -stp snpcalling2\n\t-i <input folder with mapped bam files, it will create a bamfilelist with the input unless it exist with in the inputfolder>\n\t-o <output folder to save vcf files>\n\t-rg <reference genome>\n\t-snc <number the cores to be used in parallel, recommend to use the number of Chromosomes, default 24>\n\t-chrf <a file with the list of scaffolds to use, scaffold or chromosomes names as they appear in the reference, one per row. Default N and will use the first 30 scaffold in the reference fasta file>\n\nExample:\nMELFUwgrs.pl -stp snpcalling -i ./Yuma/aligned/ -o ./Yuma/rawsnp/ -rg ./Yuma/genome/reference_genome.fasta -snc 23 -chrf chromosomefile\n\n"; exit;}
$snc //= 20;
$chrf //= "N";
if ( !-d $outputfolder ) {
    make_path $outputfolder or die "Failed to create path: $outputfolder";
}
opendir(DIR, "$inputfolder") or die "Can not open folder, $!\n" ;
our @files = readdir(DIR);
closedir (DIR);
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
chomp @scaffr;
our @cores = (1 .. (scalar @scaffr));

our $cmd="parallel --link -j $snc bcftools mpileup -Ou -r {1} -f $refgenome -b $bamlist -d 1000 -a \"FORMAT/DP,FORMAT/AD\" \'|\' bcftools call -vmO z -o $outputfolder/call.{2}.vcf.gz ::: @scaffr ::: @cores";
print "$cmd\n\n";
system($cmd);


`parallel -j $snc bcftools index $outputfolder/call.{1}.vcf.gz ::: @cores`;
`bcftools concat -a -O z --rm-dups all  $outputfolder/call.*.vcf.gz -o $outputfolder/all.vcf.gz`;
`tabix -p vcf $outputfolder/all.vcf.gz`;

}
1;
