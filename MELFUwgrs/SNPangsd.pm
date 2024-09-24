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
        elsif ($ar=~ /^-nind/){$nind=(split(/ /,$ar))[1];}
}

if (not defined ($inputfolder && $outputfolder && $refgenome)){print "\nThis script will create the ANGSD command to call and genotype SNPs for several samples in parallel, it use a likelihood approach recommended for low coverage genome sequencing. Its need the mapped bam files in a folder,and a indexed reference genome. It is recommended to use  just chromosomes or the biggest scaffolds, so you can use a file with the names of chromosomes or scaffolds to be used, default will use the any contig or scaffold with chr in the name.\n\nUsage:\nMELFUwgrs -stp snpcalling\n\t-i <input folder with mapped bam files with realigned indels>\n\t-o <output folder to save vcf files>\n\t-rg <reference genome>\n\t-snc <number the cores to be used in parallel, recommend to use the number of Chromosomes, default 20>\n\nExample:\nMELFUwgrs.pl -stp snpcalling -i ./Yuma/indelrealigned/ -o ./Yuma/rawsnp/ -rg ./Yuma/genome/reference_genome.fasta -snc 23 -chrf\n\n"; exit;}
$snc //= 20;
$chrf //= "N";
$nind //= 0.8;

opendir(DIR, "$inputfolder") or die "Can not open folder, $!\n" ;
our @files = readdir(DIR);
closedir (DIR);
if ( !-d $outputfolder ) {
    make_path $outputfolder or die "Failed to create path: $outputfolder";
}
our $LDo=$outputfolder."/LDout";
if ( !-d $LDo ) {
    make_path $LDo or die "Failed to create path: $LDo";
}
our $pruned=$outputfolder."/pruned";
if ( !-d $pruned ) {
    make_path $pruned or die "Failed to create path: $pruned";
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
our @scaffr=`grep \">\" $refgenome | head -n 30 | perl -p -e \'s/>//\' | perl -p -e \'s/\\s\.*\\n/\\n/\'`;
} else {
our @scaffr=`cat $chrf`;
}
chomp @scaffr;
$NS=scalar(@samplesnames);
$tind=int(($NS * $nind) + 0.5);
$cmd1="parallel -j $snc angsd -r {1} -b $bamlist -ref $refgenome -out $outputfolder\/{1} -GL 2 -doMajorMinor 1 -minMapQ 20 -minQ 20 -doMaf 1 -doBCF 1 -SNP_pval 1e-6 -doCounts 1 -minMaf 0.03 -doGeno 3 -doGlf 2 -uniqueonly 1 -remove_bads 1 -only_proper_pairs 0 -C 50 -baq 1 -doPost 1 -postCutoff 0.9 -geno_minDepth 3 -geno_maxDepth 1000 -nThreads 3 -minInd $tind ::: @scaffr";
##$cmd1="parallel -j $snc angsd -r {1} -b $bamlist -ref $refgenome -out $outputfolder\/{1} -GL 2 -doMajorMinor 1 -minMapQ 20 -minQ 20 -doMaf 1 -doBCF 1 -SNP_pval 1e-6 -doCounts 1 -minMaf 0.03 -doGeno 3 -doGlf 2 -uniqueonly 1 -remove_bads 1 -only_proper_pairs 0 -C 50 -baq 1 -doPost 1 -postCutoff 0.9 -geno_minDepth 3 -geno_maxDepth 1000 -nThreads 3 -minInd 17 ::: @scaffr";
print "$cmd1\n";
system ($cmd1);





use Parallel::Loops;
our %CHRNSNP =();
our $cpuCHR = scalar(@scaffr);
my $plchr = Parallel::Loops->new($cpuCHR);
$plchr->share( \%CHRNSNP);
$plchr->foreach (\@scaffr, sub{
my $chr = $_ ;
`zcat $outputfolder\/$chr\.mafs.gz | awk 'NR>1 {print \$1\"\\t\"\$2}' | gzip \>\> $outputfolder\/$chr\.pos.gz`;
`zcat $outputfolder\/$chr\.beagle.gz | cut -f 4-  |  gzip \>\> $outputfolder\/$chr\_Y.beagle.gz`;
$SNPN=`zcat $outputfolder\/$chr\.pos.gz | wc -l`;
`zcat  $outputfolder\/$chr\.mafs.gz | cut -f5 |sed 1d >>  $outputfolder\/$chr\.freq`;
`zcat $outputfolder\/$chr\.mafs.gz | awk -v OFS='\\t' 'NR>1 {print \$1,\$2,\$3,\$4}' >> $outputfolder\/$chr\_snps.list`;
`angsd sites index $outputfolder\/$chr\_snps.list`;
`angsd -r $chr -b $bamlist -ref $refgenome -out $outputfolder\/$chr -GL 2 -doGlf 3 -doMajorMinor 3 -doMAF 1 -sites $outputfolder\/$chr\_snps.list -nThreads 3`;

chomp $SNPN;
print "$SNPN\n";
$NS=scalar(@samplesnames);
`ngsLD --geno $outputfolder\/$chr\_Y.beagle.gz --pos $outputfolder\/$chr\.pos.gz --probs --n_ind $NS --n_sites $SNPN --max_kb_dist 50 --n_threads 3 --out $LDo\/$chr\.ld `;
$CHRNSNP{$chr}=$SNPN;
##For old version of ngsLD with not headers in the ld files
####`prune_graph --in $LDo\/$chr\.ld --weight-field column_7 --weight-filter \"column_3 <= 50000 && column_7 >= 0.8\" --out $LDo\/$chr\_unlinked.pos`;
##For new version of ngsLD with headers in the output ld files
`prune_graph --header --in $LDo\/$chr\.ld --weight-field \"r2\" --weight-filter \"dist <= 50000 && r2 >= 0.8\" --out $LDo\/$chr\_unlinked.pos`;
`cat $LDo\/$chr\_unlinked.pos | while read i\; do POS=\$(echo \$i | awk -F\"\:\" \'{print \$2}\')\; zcat $outputfolder\/$chr\.mafs.gz | awk -v pop=\$POS -v OFS=\"\\t\" '\$2 == pop {print \$1,\$2,\$3,\$4 ; exit}' >> $LDo\/$chr\_snps.list;  done`;
`awk -F\"\:\" \'{print \$2}\' $LDo\/$chr\_unlinked.pos | while read -r POS; do zcat $outputfolder\/$chr\.mafs | awk -v pop=\$POS -v OFS=\"\\t\" '\$2 == pop {print \$1,\$2,\$3,\$4 ; exit}' >> $LDo\/$chr\_snps.list;  done`;
`angsd sites index $LDo\/$chr\_snps.list`;
`angsd -r $chr -b $bamlist -ref $refgenome -out $pruned\/$chr -GL 2 -doGlf 2 -doMajorMinor 3 -doMAF 1 -doBCF 1 -doGeno 3 -doPost 1 -doCounts 1 -doIBS 1 -sites $LDo\/$chr\_snps.list -nThreads 3`;
`angsd -r $chr -b $bamlist -ref $refgenome -out $pruned\/$chr -GL 2 -doGlf 3 -doMajorMinor 3 -doMAF 1 -sites $LDo\/$chr\_snps.list -nThreads 3`;

});

foreach our $fold ($outputfolder,$pruned){
`zcat  $fold\/$scaffr[0]\.beagle.gz | head -n 1 | gzip >> $fold\/Somatic.beagle.gz`;
`zcat  $fold\/$scaffr[0]\.mafs.gz | head -n 1 | gzip >> $fold\/Somatic.mafs.gz`;

foreach my $chr (@scaffr){
        next if ($chr =~ m/X/ || $chr =~ m/Y/ || $chr =~ m/MT/);
        `zcat $fold\/$chr\.beagle.gz | tail -n +2 | gzip >> $fold\/Somatic.beagle.gz`;
        `zcat $fold\/$chr\.mafs.gz | tail -n +2 | gzip >> $fold\/Somatic.mafs.gz`;
        `zcat $fold\/$chr\.glf.gz | gzip >> $fold\/Somatic.glf.gz`;
}

`zcat  $fold\/Somatic.mafs.gz | cut -f6 |sed 1d >>  $fold\/Somatic.freq`;
`zcat $fold\/Somatic.mafs.gz | awk -v OFS='\\t' 'NR>1 {print \$1,\$2,\$3,\$4}' >> $fold\/Somatic_snps.list`;
`angsd sites index $fold\/Somatic_snps.list`;
}
}

1;







