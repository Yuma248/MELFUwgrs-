#!/usr/bin/perl
use File::Path qw( make_path );
while (@ARGV){
	$_=shift @ARGV;
	if ($_=~ /^-i$/){$input=shift @ARGV;}
	elsif($_=~ /^-o$/){$output=shift @ARGV;}
	elsif($_=~ /^-lnc$/){$lnc=shift @ARGV;}
	elsif($_=~ /^-cl$/){$cl=shift @ARGV;}
	elsif($_=~ /^-cn$/){$cn=shift @ARGV;}

}
if (not defined ($input)){print "\nThis pipeline calculates the depth coverage per site per sample. The required arguments are an input folder containing sam or bam files.\n\nUsage:\nExtCoverage.pl\n\t-i <Path to input folder with SAM or BAM files>\n\t-o <Path to ouput folder to save results>\n\t-lnc <Number of parallel jobs to run, default 10>\n\t-cl <File with list of scaffolds to consider; recommended: largest and closest to chromosome level>\n\t-cn <Number of scaffolds or chromosome to consider if -cl is not defined, default 24>\n\nExample:\nExtCoverage.pl -i /my/mapped/samples/ -o /my/mapped/depth/ -lnc 60 -cn 24"; exit;}
$output //= "./coverage/";
$lnc //= 10;
$cn //= 24;
if ( !-d $output) {
        make_path $output or die "Failed to create path: $output";
}
opendir(DIR, "$input") or die "Can not open folder, $!\n" ;
my @files = readdir(DIR);
closedir (DIR);
our @names=();
foreach $file (@files){next unless $file =~ /.bam$/; chomp $file; $name=$file; $name=~ s/.bam$//g; push (@nms, $name);}
if (defined ($cl)){$cn = `wc -l $cl | awk '{print \$1}'`};
$cn1 = $cn + 1;

my $cmd="parallel -j $lnc samtools coverage $input\/{1}.bam -o $output\/{1}.coverage ::: @nms";
#print "$cmd\n\n";
`parallel -j $lnc samtools coverage $input\/{1}.bam -o $output\/{1}.coverage ::: @nms`;

if (not defined ($cl)){
	`for i in \$(ls $output\/*.coverage); do samp=\$(basename "\$i"); echo -n \$samp >> $output/AverageCoverage; awk 'BEGIN{sum=0; chr=0} {if (NR>1 && NR<$cn) {sum += \$7; chr++}} END{print "\t"sum/chr}' \$i >> $output\/AverageCoverage; done`;
	`plot_coverage.R $output $cn`;
}
else {
	our $chrdir=$output."/CHR";
	if ( !-d $chrdir ) {`mkdir $chrdir`;}
	`for i in \$(ls $output\/*coverage); do samp=\$(basename "\$i"); cat $cl | while read l; do awk -v sca=\$l '\$1==sca' \$i >> $chrdir\/\$samp; done; done`;
	`for i in \$(ls $chrdir\/*.coverage); do samp=\$(basename "\$i"); echo -n \$samp >> $chrdir\/AverageCoverage; awk 'BEGIN{sum=0; chr=0} {if (NR>1) {sum += \$7; chr++}} END{print "\t"sum/chr}' \$i >> $chrdir\/AverageCoverage; done`;
	`plot_coverage.R $chrdir $cn`;
}
