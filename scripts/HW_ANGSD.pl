#!/usr/bin/perl
while (@ARGV){
	$_=shift @ARGV;
	if ($_=~ /^-i$/){$input=shift @ARGV;}
	elsif($_=~ /^-o$/){$output=shift @ARGV;}
	elsif($_=~ /^-b$/){$bamlist=shift @ARGV;}
	elsif($_=~ /^-p$/){$popmap=shift @ARGV;}
	elsif($_=~ /^-hwe$/){$hwe=shift @ARGV;}
	elsif($_=~ /^-rg$/){$rg=shift @ARGV;}
	elsif($_=~ /^-lnc$/){$lnc=shift @ARGV;}
	elsif($_=~ /^-cl$/){$chrlist=shift @ARGV;}
	elsif($_=~ /^-pc$/){$pcutoff=shift @ARGV;}
}
if(not defined ($input && $bamlist && $popmap && $rg && $chrlist )){print "This script use angsd to filtere SNPs that violate HWE in a proportion of the populations. Requires a input folder with the snp_list per chromosome, a bamlits file with all the bam files used,a popmap with the samples and populations defined and a chomosome list.\n\nUsages:\nHW_ANGSD.pl\n\t-i <Path to input folder with snp_list>\n\t-b <Path to file with original bam files list>\n\t-rg <Path to reference genome>\n\t-p <Path to popmap file>\n\t-cl <Path to file with the chromosomes or scaffolds to analysis>\nOptional:\n\t-o <Path to output folder, default ./HWESNP>\n\t-hwe <Pvalue threshold, defaul 0.05>\n\t-pc <Population cutoff, proportion of pop in which loci can be out of HWE, defaul 0.1>\n\t-lnc <Runs in parallel, default 24>\n\nExample:\n\nHW_ANGSD.pl -i /LDout/ -o /HWESNP -b /path/bamlist -p /path/popmap -rg /path/referencegenome -cl /paht/chromosomeslist -hwe 0.01 -pc 0.2 -lnc 24\n\n";exit;}
$output //="./HWESNP";
$lnc //= 24;
$hwe //= 0.05;
$pcutoff //=0.1;
use Parallel::Loops;

if ( !-d $output ) {
    make_path $output or die "Failed to create path: $output";
}


open my $fh, '<', $chrlist or die "Can't open $chrlist: $!";
my @scaffr = <$fh>;        # Slurp all lines (keeps newlines)
close $fh;
chomp(@scaffr);


open(POP, "<", $popmap) or die $!;
my %pops;
while(<POP>) {
	next if $_ =~ /^\s/;
	chomp;
	my ($sample, $pop) = split;
	$pops{$pop} = [] unless $pops{$pop};
	push @{$pops{$pop}}, $sample;
}
close POP;


my %exclude_count;
foreach my $pop (sort keys %pops) {
	my $indfile = $output . '/' . $pop . '.bamlist';
	open(INDO, ">", $indfile) or die $!;
	foreach my $ind (@{$pops{$pop}}) {
		my $bamfile=`grep "$ind\_" $bamlist`;
		print INDO $bamfile;
	}
	close INDO;

	#my $cmd="parallel -j $lnc angsd -r {1} -b $indfile -ref $rg -out $output\/{1}_$pop -GL 2 -doHWE 1 -domajorminor 1 -maxHWEpval $hwe -sites $input\/{1}\_snps.list -nThreads 3 ::: @scaffr";
	#print "$cmd\n";
	my $plchr = Parallel::Loops->new($lnc);
	$plchr->share(\%exclude_count);
	$plchr->foreach(\@scaffr, sub{
		my $chr = $_;
		my $cmd ="angsd -r $chr -b $indfile -ref $rg -out $output\/$chr\_$pop -GL 2 -doHWE 1 -domajorminor 1 -maxHWEpval $hwe -sites $input\/$chr\_snps.list -nThreads 3";
		system($cmd) == 0 or die "Comman $cmd failed: $?"; 
		my $HWF = "$output\/$chr\_$pop\.hwe.gz";
		open(HWEI, "-|","zcat $HWF" ) or die "Can not zcat $HWF: $!";
		<HWEI>;
		while(<HWEI>) {
			last if $_ =~ /^\s/;
			chomp;
			my ($locus, $pos) = split;
			$exclude_count{"$locus-$pos"}++;
		}
		close HWEI;
		
	});


}

my $plchr2 = Parallel::Loops->new($lnc);
$plchr2->foreach(\@scaffr, sub{
	my $chr2 = $_;
	my $HWEOF = "$output\/$chr2\_snps.list";
	open (HWEO, ">", $HWEOF) or die "Can not open $HWEOF: $!";
	my $SNPOF= "$input\/$chr2\_snps.list";
	open(SNPO, "<", $SNPOF ) or die "Can not open $SNPOF: $!";
	while(<SNPO>) {
		last if $_ =~ /^\s/;
		chomp;
		my ($locus, $pos) = split;
		my $snp = "$locus\-$pos";
		my $count = $exclude_count{$snp} || 0;
		if ($count / scalar(keys %pops) < $pcutoff) { print HWEO $_ , "\n";}
	}
	close SNPO;
	close HWEO;
	`angsd sites index $HWEOF`;
	`angsd -r $chr2 -b $bamlist -ref $rg -out $output\/$chr2 -GL 2 -doGlf 2 -doMajorMinor 3 -doMAF 1 -doBCF 1 -doGeno 3 -doPost 1 -doCounts 1 -doIBS 1 -sites $output/$chr2\_snps.list -nThreads 3`;
	`angsd -r $chr2 -b $bamlist -ref $rg -out $output\/$chr2 -GL 2 -doGlf 3 -doMajorMinor 3 -doMAF 1 -sites $output\/$chr2\_snps.list -nThreads 3`;
	
});

`zcat  $output\/$scaffr[0]\.beagle.gz | head -n 1 | gzip >> $output\/Somatic.beagle.gz`;
`zcat  $output\/$scaffr[0]\.mafs.gz | head -n 1 | gzip >> $output\/Somatic.mafs.gz`;
foreach my $chr (@scaffr){
        `zcat $output\/$chr\.beagle.gz | tail -n +2 | gzip >> $output\/Somatic.beagle.gz`;
        `zcat $output\/$chr\.mafs.gz | tail -n +2 | gzip >> $output\/Somatic.mafs.gz`;
        `zcat $output\/$chr\.glf.gz | gzip >> $output\/Somatic.glf.gz`;
}
`zcat  $output\/Somatic.mafs.gz | cut -f6 |sed 1d >>  $output\/Somatic.freq`;
`zcat $output\/Somatic.mafs.gz | awk -v OFS='\\t' 'NR>1 {print \$1,\$2,\$3,\$4}' >> $output\/Somatic_snps.list`;
`angsd sites index $output\/Somatic_snps.list`;
