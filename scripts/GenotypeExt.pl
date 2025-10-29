#!/usr/bin/perl -w
while (@ARGV){
	$_=shift @ARGV;
	if ($_=~ /^-i$/){$inputfolder=shift @ARGV;}
	elsif($_=~ /^-o$/){$outputfolder=shift @ARGV;}
	elsif ($_=~ /^-snc$/){$snc=shift @ARGV;}
	elsif ($_=~ /^-rg$/){$reference=shift @ARGV;}
	elsif ($_=~ /^-bam$/){$bamlist=shift @ARGV;}
	elsif ($_=~ /^-chrf$/){$chrf=shift @ARGV;}
}


if (not defined ($inputfolder && $outputfolder && $refgenome && $bamlist)){print "\nThis script will create the ANGSD command to genotype SNPs for extra samples when you already call SNP from another set of samples, it run in parallel, it use a likelihood approach recommended for low coverage genome sequencing. Its need the bam files list, a indexed reference genome, and the LDout folder with the _snps.list files, or the prune folder with the Somatic_snp.list file. It is recommended to use  just chromosomes or the biggest scaffolds, so you can use a file with the names of chromosomes or scaffolds to be used, default will use the any contig or scaffold with chr in the name.\n\nUsage:\nGenotypeExt.pl\n\t-i <input folder with _snp.list files>\n\t-o <output folder to save all the produced files>\n\t-rg <reference genome>\n\t-snc <number the cores to be used in parallel, recommend to use the number of Chromosomes, default 24>\n\nExample:\nGenotypeExt.pl -i ./Yuma/LDout/ -o ./Yuma/ALLSAMPLES/ -rg ./Yuma/genome/reference_genome.fasta -snc 24 -chrf my chrfile -bam mybamlist\n\n"; exit;}
$snc //= 24;
$chrf //= "N";
#$nind //= 0.8;
$sxchr //= "SEXCHROM";
$sychr //= "SEYCHROM";
$mtchr //= "MITCHROM";

$LDo //= $inputfolder;
$fold //= $outputfolder;
if ( !-d $outputfolder ) {
    make_path $outputfolder or die "Failed to create path: $outputfolder";
}
if (-e $bamlist){
    @samplesnames=`awk -F'/' '{print $NF}' $bamlist`;
   chomp(@samplesnames);
}
if ($chrf eq "N"){
our @scaffr=`grep \">\" $refgenome | head -n 30 | perl -p -e \'s/>//\' | perl -p -e \'s/\\s\.*\\n/\\n/\'`;
} else {
our @scaffr=`cat $chrf`;
}
chomp @scaffr;
if ($LDo =~ /LDout/){
use Parallel::Loops;
our %CHRNSNP =();
our $cpuCHR = scalar(@scaffr);
my $plchr = Parallel::Loops->new($cpuCHR);
$plchr->share( \%CHRNSNP);
$plchr->foreach (\@scaffr, sub{
my $chr = $_ ;
`angsd -r $chr -b $bamlist -ref $refgenome -out $fold\/$chr -GL 2 -doGlf 2 -doMajorMinor 3 -doMAF 1 -doBCF 1 -doGeno 3 -doPost 1 -doCounts 1 -doIBS 1 -sites $LDo\/$chr\_snps.list -nThreads 3`;
`angsd -r $chr -b $bamlist -ref $refgenome -out $fold\/$chr -GL 2 -doGlf 3 -doMajorMinor 3 -doMAF 1 -sites $LDo\/$chr\_snps.list -nThreads 3`;
});

`zcat  $fold\/$scaffr[0]\.beagle.gz | head -n 1 | gzip >> $fold\/Somatic.beagle.gz`;
`zcat  $fold\/$scaffr[0]\.mafs.gz | head -n 1 | gzip >> $fold\/Somatic.mafs.gz`;

foreach my $chr (@scaffr){
#        next if ($chr =~ m/X/ || $chr =~ m/Y/ || $chr =~ m/MT/);
    next if ($chr eq $sychr || $chr eq $sxchr || $chr eq $mtchr);
    `zcat $fold\/$chr\.beagle.gz | tail -n +2 | gzip >> $fold\/Somatic.beagle.gz`;
    `zcat $fold\/$chr\.mafs.gz | tail -n +2 | gzip >> $fold\/Somatic.mafs.gz`;
    `zcat $fold\/$chr\.glf.gz | gzip >> $fold\/Somatic.glf.gz`;
}

`zcat  $fold\/Somatic.mafs.gz | cut -f6 |sed 1d >>  $fold\/Somatic.freq`;
`zcat $fold\/Somatic.mafs.gz | awk -v OFS='\\t' 'NR>1 {print \$1,\$2,\$3,\$4}' >> $fold\/Somatic_snps.list`;
`angsd sites index $fold\/Somatic_snps.list`;
}
elsif($LDo =~ /pruned/){
use Parallel::Loops;
our %CHRNSNP =();
our $cpuCHR = scalar(@scaffr);
my $plchr = Parallel::Loops->new($cpuCHR);
$plchr->share( \%CHRNSNP);
$plchr->foreach (\@scaffr, sub{
my $chr = $_ ;
`angsd -r $chr -b $bamlist -ref $refgenome -out $fold\/$chr -GL 2 -doGlf 2 -doMajorMinor 3 -doMAF 1 -doBCF 1 -doGeno 3 -doPost 1 -doCounts 1 -doIBS 1 -sites $LDo\/Somatic_snps.list -nThreads 3`;
`angsd -r $chr -b $bamlist -ref $refgenome -out $fold\/$chr -GL 2 -doGlf 3 -doMajorMinor 3 -doMAF 1 -sites $LDo\/Somatic_snps.list -nThreads 3`;
});

`zcat  $fold\/$scaffr[0]\.beagle.gz | head -n 1 | gzip >> $fold\/Somatic.beagle.gz`;
`zcat  $fold\/$scaffr[0]\.mafs.gz | head -n 1 | gzip >> $fold\/Somatic.mafs.gz`;

foreach my $chr (@scaffr){
#        next if ($chr =~ m/X/ || $chr =~ m/Y/ || $chr =~ m/MT/);
    next if ($chr eq $sychr || $chr eq $sxchr || $chr eq $mtchr);
    `zcat $fold\/$chr\.beagle.gz | tail -n +2 | gzip >> $fold\/Somatic.beagle.gz`;
    `zcat $fold\/$chr\.mafs.gz | tail -n +2 | gzip >> $fold\/Somatic.mafs.gz`;
    `zcat $fold\/$chr\.glf.gz | gzip >> $fold\/Somatic.glf.gz`;
}

`zcat  $fold\/Somatic.mafs.gz | cut -f6 |sed 1d >>  $fold\/Somatic.freq`;
`zcat $fold\/Somatic.mafs.gz | awk -v OFS='\\t' 'NR>1 {print \$1,\$2,\$3,\$4}' >> $fold\/Somatic_snps.list`;
`angsd sites index $fold\/Somatic_snps.list`;
}
