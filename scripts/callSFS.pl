#!/usr/bin/perl -w

use Parallel::Loops;

while (@ARGV){
        $_=shift @ARGV;
        if ($_=~ /^-pm$/){$pm=shift @ARGV;}
        elsif($_=~ /^-o$/){$output=shift @ARGV;}
        elsif ($_=~ /^-i$/){$input=shift @ARGV;}
        elsif ($_=~ /^-rg$/){$reference=shift @ARGV;}
        elsif ($_=~ /^-lnc$/){$lnc=shift @ARGV;}
	elsif ($_=~ /^-s$/){$lsnp=shift @ARGV;}
	elsif ($_=~ /^-snc$/){$snc=shift @ARGV;}
}
if (not defined ($pm && $output && $input && $reference)){print "\nThis script uses angsd to calculate SFS per indiviudal and population, and then calculate He, and Fst base on these SFS. The required arguments are popmap that assigns individual to each pop, an input folder with the bam files of all the samples, and a output folder to save all the outputs.\nUsage:\ncalSFS\n\t-i <inputfolder with all bam files to used>\n\t -o <output folder to save all the resutls>\n\t -pm <popmap, tab delimited files with sample names in the first column and corresponding pop in the second column>\n\t -rg <path to reference genome used to map the bam files>\nOptional:\n\t-s <liste of sites to be considere, recomended to use pruned SNPs to avoid LD effects, this list should be in the site format from ANGSD>\n\t-lnc <number of runs in parallel, it works per sample pop and pairs of pop, default 10>\n\t-snc <number of CPU to use in pairwise Fst calculations, default 4>\n\nExample:\ncalSFS -i /Yuma/mapped/ -o /Yuma/SFSout/ -pm /Yuma/popmap -rg /Yuma/genome/REFERENCE.fasta -s /Yuma/pruned/Somatic_prunedSNPs.list -snc 4 -lnc 25.\n\n";exit;}
if (not defined ($lnc)){$lnc = 10};
if (not defined ($snc)) {$snc = 4};

our @maps = `cat $pm`;
our @pops = `awk '{print \$2}' $pm | sort | uniq`;
our @Abam = ();
chomp @pops;


my $plpop = Parallel::Loops->new($lnc);
$plpop->foreach (\@pops, sub{
my $pop = $_ ;
	my $poplist=$output."/".$pop."_list";
	if (-e $poplist) {
		our @Pbam=`cat $poplist`;
		push (@Abam, @Pbam);
	}else{
		our @Pbam=();
		my @samples = map { (split /\t/)[0] } grep { /$pop$/ } @maps;
		for $sample (@samples){
			my $bam = $input.$sample."_repeatmasked.bam";
			`echo $bam >> $poplist`;
			push(@Pbam, $bam);
			push(@Abam, $bam);
		}	
	}
	if (defined ($lsnp)){
		our $cmd = "angsd -b $poplist -anc $reference -out $output\/$pop -doCounts 1 -GL 1 -doSaf 1 -sites $lsnp";
	}else{
		our  $cmd = "angsd -b $poplist -anc $reference -out $output\/$pop -doCounts 1 -GL 1 -doSaf 1";
	}
	system($cmd); 	
	our $cmd2 = "realSFS $output\/$pop\.saf.idx -fold 1 >> $output\/$pop\.folded.sfs";
	system ($cmd2);
	our $cmd5 ="realSFS saf2theta $output\/$pop\.saf.idx -sfs $output\/$pop\.folded.sfs -outname $output\/$pop\_Theta";
	system ($cmd5);
	our $cmd6 = "thetaStat do_stat $output\/$pop\_Theta.thetas.idx -win 10000 -step 10000 -outnames $output\/$pop\.thetas.windows.gz";
	system ($cmd6);
	our $cmd7 = "awk -F '\t' '\$14>0 && NR>1  {print \$2, \$3, \$4 / \$14, \$5 / \$14 }' $output\/$pop\.thetas.windows.gz.pestPG  >> $output\/$pop\.corrected.thetas.pi";
	system ($cmd7);
});

my $plpop3 = Parallel::Loops->new($lnc);
$plpop3->foreach (\@pops, sub{
my $pop = $_ ;
	my $plpop2 = Parallel::Loops->new($lnc);
	$plpop2->foreach (\@pops, sub{
		my $pop2 = $_ ;
		unless ($pop2 eq $pop){
			our $cmd3 = "realSFS $output\/$pop\.saf.idx $output\/$pop2\.saf.idx >> $output\/$pop\.$pop2\.sfs";
			system ($cmd3);
			our $cmd4 = "realSFS fst index $output\/$pop\.saf.idx $output\/$pop2\.saf.idx -sfs $output\/$pop\.$pop2\.sfs -fstout $output\/$pop\_$pop2\_fst -whichFst 1 -P $snc";
			system ($cmd4);
			our $cmd42 = "realSFS fst stats $output\/$pop\_$pop2\_fst.fst.idx >> $output\/$pop\_$pop2\_FST";
			system ($cmd42);
		}
	});	
});

my $cmd8 = "for POP in @pops ; do echo -n \${POP}\$'\\t' >> $output\/THETA_All; awk -v OFS='\\t' '{ sum += \$3; sum2 += \$4; n++ } END { if (n > 0) print sum / n, sum2 /n;}' $output\/\${POP}.corrected.thetas.pi >> $output\/THETA_All; done";
system ($cmd8);
my $cmd9 = "for POP in @pops ;  do echo -n \${POP}\$'\\t' >> $output\/HeAll;  awk -v OFS='\\t' '{for(i=1;i<=NF;i++) t+=\$i; print \$2/t; t=0}' $output\/\${POP}.folded.sfs >> $output\/HeAll; done";
system ($cmd9);
my $cmd10 = "for POP in @pops; do for POP2 in @pops; do if [[ \$POP != \$POP2 ]]; then echo -n \${POP}\$'\\t'\${POP2}\$'\\t' >> $output\/All_FSTs; cat $output\/\${POP}_\${POP2}_FST >> $output\/All_FSTs; else continue; fi;  done; done"; 
system ($cmd10);
my $tnc = $lnc * $snc;
our @sams = `awk '{print \$1}' $pm | sort | uniq`;
chomp @sams;
my $plsam = Parallel::Loops->new($tnc);
$plsam->foreach (\@sams, sub{
	my $sam = $_ ;
	my $bam  = $input.$sam."_repeatmasked.bam";
		if (defined ($lsnp)){
		our $cmd11 = "angsd -i $bam -anc $reference -out $output\/$sam -doCounts 1 -GL 1 -doSaf 1 -sites $lsnp";
		}else{
		our  $cmd11 = "angsd -i $bam -anc $reference -out $output\/$sam -doCounts 1 -GL 1 -doSaf 1";
		}
	system($cmd11);
	our $cmd12 = "realSFS $output\/$sam\.saf.idx -fold 1 >> $output\/$sam\.folded.sfs";
	system ($cmd12);
});
my $cmd13 = "for SAM in @sams ; do echo -n \${SAM}\$'\\t' >> $output\/HeAllS;  awk -v OFS='\\t' '{for(i=1;i<=NF;i++) t+=\$i; print \$2/t; t=0}' $output\/\${SAM}.folded.sfs >> $output\/HeAllS; done";
system ($cmd13);




