#!/usr/bin/perl -w
package NGSadmix;

my $LEVEL = 1;
sub NGSadmixP{
my @arg = @_;
foreach $ar (@arg){
	if  ($ar =~ /^-i/){our $input = (split(/ /,$ar))[1];}
        elsif($ar=~ /^-o/){our $output = (split(/ /, $ar))[1];}
        elsif ($ar=~ /^-nc/){$nc=(split(/ /, $ar))[1];}
        elsif ($ar=~ /^-snc/){$snc=(split(/ /, $ar))[1];}
	elsif ($ar=~ /^mK$/){$mK=(split(/ /, $ar))[1];}
	elsif ($ar=~ /^-rep$/){$rep=(split(/ /, $ar))[1];}
 	elsif ($ar=~ /^-ind$/){$ind=(split(/ /, $ar))[1];}
	elsif ($ar=~ /^-maf$/){$maf=(split(/ /, $ar))[1];}
}

if (not defined ($input )){print "\nUsage:\nNGSadmixP -i <beagle format input file>\nOptional:\n\t-o <output folder to save resutls, default ./NGSadmixR>\n\t-ind <minimum number of informative individuals, default 0>\n\t-mK <maximum K, number of populations assumed, default 5>\n\t-nc <number cores to run in parallel, default 2 >\n\t-snc <number of cores per each run, default 30>\n\t-rep <number replicates, default 1>\n\t-maf <minimum minor allele frequency, default 0.01>\n\nFor example:\n NGSadmixP -i yuma.str -o ./NGSadmixR/ -mK 10 -rep 5 -nc 3 -snc 30 -maf 0.05 \n\n"; exit;} 
if (not defined $rep){$rep = 1};
if (not defined $mK){$mK = 5};
if (not defined $nc){$nc = 2};
if (not defined $snc){$snc = 30};
if (not defined $maf){$maf = 0.01};
if (not defined $ind){$ind = 0};
if (not defined $output){$output = "./NGSadmixR/"};
if (! -d $output){
	`mkdir $output`;
}
our @Ks = (1..$mK);
our @Rs = (1..$rep);

`parallel -j $nc NGSadmix -likes $input -minMaf $maf -K {1} -o $output\/NGSadmix_$ind\_{1}_{2} -minInd $ind -P $snc ::: @Ks ::: @Rs `;

}
1;
