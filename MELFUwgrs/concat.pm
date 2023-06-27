package concat;

my $LEVEL = 1;
sub coca{
my @arg = @_;
foreach $ar (@arg){
	if ($ar =~ /^-i/){our $inputfolder = (split(/ /,$ar))[1];}
	elsif($ar=~ /^-o/){our $outputfolder = (split(/ /, $ar))[1];}
	elsif ($ar=~ /^-snc/){$snc=(split(/ /, $ar))[1];}
	elsif ($ar=~ /^-t/) {$type = (split(/ /, $ar))[1];}
	elsif ($ar=~ /^-exf/) {$exf = (split(/ /, $ar))[1];}
}
use File::Path qw( make_path );
if (not defined ($inputfolder && $outputfolder)){print "\nThis script will concatenate fastq files of raw or trimmed reads of several samples in parallel. It requires fastq files of all samples stored in the same folder\n\nUsage:\nConcat.pl\n\t-i <path to the folder with the input fastq files>\n\t-o <path to the output folder>\n Optional:\n\t-snc <number of runs in parallel, default 10>\n\t-t <method used for trimming. Trimmomatic TR, AdapterRemoval AR or None NO if are raw sequences, default AR>\n\t-exf <this will tell the script how to extract the name of each sample, and should include all extra information at the end of the file names that is not related to the sample name, default P1_L001_. >\n\nFor example:\nconcat.pl -i /home/Yumafan/demultiplex/trimmed/ -o /home/Yumafan/concatenated-snc 10 -t AR -exf P1_L001_,\n\n"; exit;}


if ( !-d $outputfolder ) {
	make_path $outputfolder or die "Failed to create path: $outputfolder";
}
our $tmpdir=$outputfolder."/tmp";
if ( !-d $tmpdir ) {
        `mkdir $tmpdir`;
}
if (not defined ($snc)){$snc =10;}

if ($type eq "AR"){
        our $code =$exf.".collapsed.gz";
        my @names = `ls $inputfolder\/*$code`;
        foreach $name (@names) {chomp $name; $name=~ s/$inputfolder\///g; $name=~ s/$code//g; push (@nms, $name);}
        foreach $ef ( "singleton.truncated", "collapsed", "collapsed.truncated", "pair1.truncated", "pair2.truncated"){
                my $cmd = "parallel -j $snc --results $tmpdir --tmp $tmpdir zcat $inputfolder\/{1}*\.$ef\.gz | gzip \'>>\' $outputfolder/{1}\.$ef\.gz ::: @nms";
                print "$cmd\n";
                `parallel -j $snc --results $tmpdir --tmp $tmpdir zcat $inputfolder\/{1}*\.$ef\.gz \'|\' gzip \'>\'\'>\' $outputfolder/{1}\.$ef\.gz ::: @nms`
                }
}
elsif ($type eq "TR"){
        my $code = $exf.".gz";
        my @names = `ls $inputfolder\/\*$code`;
        foreach $name (@names) {chomp $name; $name=~ s/$inputfolder\///g; $name=~ s/$code//g; push (@nms, $name);}
        my $cmd = "parallel -j $snc --results $tmpdir --tmp $tmpdir zcat $inputfolder\/{1}*\.gz | gzip \'>>\' $outputfolder/{1}\.gz ::: @nms";
        print "$cmd\n";
        `parallel -j $snc --results $tmpdir --tmp $tmpdir zcat $inputfolder\/{1}*\.gz \'|\' gzip \'>\'\'>\' $outputfolder/{1}\.gz ::: @nms`
}
elsif ($type eq "NO"){
        my $code = $exf.".gz";
        my @names = `ls $inputfolder\/\*$code`;
        foreach $name (@names) {chomp $name; $name=~ s/$inputfolder\///g; $name=~ s/$code//g; push (@nms, $name);}
	foreach $ef ( "1", "2"){
                my $cmd = "parallel -j $snc --results $tmpdir --tmp $tmpdir zcat $inputfolder\/{1}_$ef\.fq.gz | gzip \'>>\' $outputfolder/{1}\.$ef\fq.gz ::: @nms";
                print "$cmd\n";
                `parallel -j $snc --results $tmpdir --tmp $tmpdir zcat $inputfolder\/{1}_$ef\.fq.gz \'|\' gzip \'>\'\'>\' $outputfolder/{1}\.$ef\.fq.gz ::: @nms`
                }
else {print "You have to select a correct type of file, check the instrcutions\n"}


} 

1;
