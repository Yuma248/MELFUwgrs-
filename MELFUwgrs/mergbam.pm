package mergbam;

my $LEVEL = 1;
sub merg{
my @arg = @_;
foreach $ar (@arg){
        if ($ar =~ /^-i /){our $inputfolder = (split(/ /,$ar))[1];}
        elsif($ar=~ /^-o /){our $outputfolder = (split(/ /, $ar))[1];}
        elsif ($ar=~ /^-snc /){$snc=(split(/ /, $ar))[1];}
        elsif ($ar=~ /^-exf /) {$exf = (split(/ /, $ar))[1];}
        elsif ($ar=~ /^-ind /) {$ind = (split(/ /, $ar))[1];}
}
use File::Path qw( make_path );
if (not defined ($inputfolder && $outputfolder)){print "\nThis script will merge bam files of different lanes or runs of the same sample, it is important that you check the runs have similar quality and have the same category of read (see samtools flagstat). The script can work wiht several samples in parallel. It requires bam files of all samples stored in the same folder\n\nUsage:\nSNPcallPipe.pl -stp merge\n\t-i <path to the folder with the input bam files>\n\t-o <path to the output folder>\n Optional:\n\t-snc <number of runs in parallel, default 10>\n\t-exf <this will tell the script how to extract the name of each sample, and should include all extra information at the end of the file names that is not related to the sample name, default _S*P1_L001_>\n\t-ind <if you need to samtools index the input files, default yes>\n\nFor example:\nSNPcallconcat.pl -i /home/Yumafan/dedupout/ -o /home/Yumafan/mergedbam -snc 10 -exf S*P1_L001__markdup -ind yes\n\n"; exit;}
if ( !-d $outputfolder ) {
        make_path $outputfolder or die "Failed to create path: $outputfolder";
}
our $tmpdir=$outputfolder."/tmp";
if ( !-d $tmpdir ) {
        `mkdir $tmpdir`;
}
if (not defined ($snc)){$snc =10;}
if (not defined ($exf)){$exf ="_S\*P1_L001__markdup";}
our  $code = $exf.".bam";
our @names = `ls $inputfolder\/\*$code`;
$code =~ s/\*/\\d\{1,2\}/;
foreach $name (@names) {
        chomp $name;
        $name=~ s/$inputfolder\///;
        $name=~ s/$code//;
        push (@nms, $name);
}
our @uninms = grep { ! $seen{$_} ++ } @nms;
use Parallel::Loops;
my $plsam = Parallel::Loops->new($snc);
$plsam->foreach (\@uninms, sub{
my $unms = $_ ;
my @listf = `ls $inputfolder\/$unms\*\.bam`;
chomp @listf;
my $cmd2="samtools merge $outputfolder\/$unms\_mergeA.bam ";
my $cmd3="samtools index -c  $outputfolder\/$unms\_mergeA.bam";
foreach $file (@listf){
        if ($ind eq "yes"){
                my $cmd = "samtools index -c $file";
                system($cmd);
        }
        $cmd2=$cmd2.$file." ";
}
system($cmd2);
system($cmd3);

});

                                                }
1;


