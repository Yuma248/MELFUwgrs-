#!/usr/bin/perl -w

use Parallel::Loops;

while (@ARGV){
$_=shift @ARGV;
if ($_=~ /^-if$/){$input=shift @ARGV;}
elsif ($_=~ /^-pm$/){$popmap=shift @ARGV;}
elsif ($_=~ /^-of$/){$output=shift @ARGV;}
elsif ($_=~ /^-nc$/){$nc=shift @ARGV;}
}
if (not defined ($input && $popmap)){print "\nThis script calcualate pairwise Fst and population specific Fst first per SNP  and then create a table with global values, it uses vcftools for the calculations.\n\nUsage:\nFstVCFtools.pl\n\t-if <path to a vcf input file>\n\t-pm <path to a popmap file, tabdelimted two columns, <ID> <pop>>\nOptional:\n\t-of <path to output folder, deafult ./outputFst>\n\t-nc <number of cores to run in parallel, deaful 5>\n\nFor example:\n\nFst_VCFtools.pl -if /ANGSD/pruned -of /SomaticFST/ -pm /lineages/popmap -nc 58\n\n"; exit;}
$output //= "./outputFst/";
mkdir $output unless -d $output;
$nc //= 5;
my @mypops;
my @pairpops;
my %pops;
my @allsamples;
open(my $POP, "<", $popmap) or die $!;
while(<$POP>) {
	next if $_ =~ /^\s/;
	chomp;
	my ($sample, $pop) = split;
	push @{ $pops{$pop} //= [] }, $sample;
	push @allsamples , $sample;
}
close POP;

foreach my $pop (sort keys %pops) {
	push @mypops, $pop;
	my $indfile = "${output}/${pop}.inds";
	my $Nindfile = "${output}/${pop}_INV.inds";
	my %indh;
	$indh{$_}=1 for @{$pops{$pop}};
	my @inv_pop = grep { !$indh{$_} } @allsamples;
	next if -e $indfile;
	open(my $INDO, ">", $indfile) or die $!;
	print $INDO "$_\n" for @{$pops{$pop}};
	close INDO;
	open (my $ININV, ">", $Nindfile) or die $!;
	print $ININV "$_\n" for @inv_pop;
	close $ININV;
}
for my $ipop (0 .. $#mypops - 1){
	for my $jpop ($ipop + 1 .. $#mypops){
		push @pairpops, "${mypops[$ipop]}_${mypops[$jpop]}";
	}
}

my $plraw = Parallel::Loops->new($nc);
$plraw->foreach (\@pairpops, sub{
	my $pair= $_ ;
	my ($fpop, $spop) = split /_/, $pair;
	`vcftools --vcf $input --weir-fst-pop ${output}/${fpop}.inds --weir-fst-pop ${output}/${spop}.inds --out ${output}/${fpop}_${spop}_FST`;


});

my $fstfile = "${output}/pairwise.fst";
open (my $INFST, ">", $fstfile) or die $!;
foreach my $pairk (@pairpops) {
	my ($fpop, $spop) = split /_/, $pairk;
	my $wfst = `grep "weight"  ${output}/${pairk}_FST.log | awk '{print \$NF}'`;
	chomp $wfst;
	my $ufst = `grep "mean"  ${output}/${pairk}_FST.log | awk '{print \$NF}'`;
	chomp $ufst;
	print $INFST "${fpop}\t${spop}\t${ufst}\t${wfst}\n";
}
close $INFST;

my $plraw2 = Parallel::Loops->new($nc);
$plraw2->foreach (\@mypops, sub{
	my $pop = $_ ;
	`vcftools --vcf $input --weir-fst-pop ${output}/${pop}.inds --weir-fst-pop ${output}/${pop}_INV.inds --out ${output}/${pop}_specific_FST`;
});

my $Sfstfile = "${output}/specific.fst";
open (my $INSFST, ">", $Sfstfile) or die $!;
foreach my $pop (@mypops) {
	my $wfst = `grep "weight"  ${output}/${pop}_specific_FST.log | awk '{print \$NF}'`;
	chomp $wfst;
	my $ufst = `grep "mean"  ${output}/${pop}_specific_FST.log | awk '{print \$NF}'`;
	chomp $ufst;
	print $INSFST "${pop}\t${ufst}\t${wfst}\n";
}


