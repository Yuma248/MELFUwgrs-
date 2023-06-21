
package indelrealigment;

my $LEVEL = 1;
sub indelreal {

use File::Path qw( make_path );
my @arg = @_;
foreach $ar (@arg){
        if ($ar=~ /^-i /){our $inputfolder= (split(/ /,$ar))[1];}
        elsif ($ar=~ /^-o /){our $outputfolder= (split(/ /,$ar))[1];}
        elsif ($ar=~ /^-ind /){our $ind=(split(/ /,$ar))[1];}
        elsif ($ar=~ /^-ncp /){our $ncp=(split(/ /,$ar))[1];}
	elsif ($ar=~ /^-rg /){$refgenome=(split(/ /,$ar))[1];}
	elsif ($ar=~ /^-lnc /){$lnc=(split(/ /,$ar))[1];}
	
}

if (not defined ($refgenome && $inputfolder && $outputfolder)){print "\nThis script will create a list of indel per sample, and it will realign reads around indels for a more efficient variants calling. It uses GATK in parallel to process several samples at the same time.\n\nUsage:\nIndelRealigment\n\t-i <input folder with sorted bam files>\n\t-o <output folder to save the realigments>\n\t-rg <path to the reference genome>\nOptional:\n\t-ncp <number cores to run in parallel, default 4>\n\t-ind <\"yes\" if you need to index your bam files, or \"no\" if they are already indexed, default \"yes\">\n\nFor example: IndelRealigment -i ./Yuma/sam2bamout/ -o ./Yuma/indelout/ -rg ./Yuma/genomes/reference.fasta -ncp 8 -ind yes\n\n"; exit;}
opendir(DIR, "$inputfolder") or die "Can not open inputfolder $inputfolder, $!\n" ;
my @files = readdir(DIR);
closedir (DIR);
if ( !-d $outputfolder ) {
    make_path $outputfolder or die "Failed to create path: $outputfolder, $!\n";
}
if (not defined ($ncp)){$ncp=8;}
if (not defined ($ind)){$ind="yes";}
our @samplesnames=();
our @pops=();
our @term=();
#foreach my $file (@files){next unless ($file =~ /\.bam$/ and $file =~ /^sort_/); my @fileinf =split (/[\_\.]/, $file); if ($fileinf[1] =~ m/^pop/){push @samplesnames,$fileinf[2];push @pops, $fileinf[1];}else {push @samplesnames, $fileinf[1];}}
foreach my $file (@files){next unless ($file =~ /\.bam$/ ); my @fileinf =split (/[\_\.]/, $file); if ($fileinf[0] =~ m/^pop/){push @samplesnames,$fileinf[1];push @pops, $fileinf[0];push @term, $fileinf[2]; }else {push @samplesnames, $fileinf[0]; push @term, $fileinf[1];}}
#foreach my $file (@files){ next unless ($file =~ /\.bam$/ and $file =~ /markdup/); my @fileinf = split (/[\_\.]/, $file); push @samplesnames, $fileinf[1];}
#foreach my $name (@samplesnames){my $forward = "$inputfolder\/$name\.trim\.1\.fq\.gz"; my $reverse = "$inputfolder\/$name\.trim\.2\.fq\.gz"; my $output = " $outputfolder\/$name\.sam"; print "$name\t$forward\t$reverse\t$output\n";}


if (@pops){
$cmd2="parallel --xapply -j $lnc --results ./logfilesIR --noswap samtools index $inputfolder\/{2}_{1}_{3}\.bam ::: @samplesnames ::: @pops ::: @term";
$cmd1="parallel --xapply -j $ncp --results ./logfilesIR --noswap java -jar ./SNPcallPipe/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $refgenome -I $inputfolder\/{2}\_{1}_{3}\.bam -o $outputfolder\/{2}\_{1}\_intervals\.list ::: @samplesnames ::: @pops ::: @term";
$cmd="parallel --xapply -j $ncp --results ./logfilesIR --noswap java -jar ./SNPcallPipe/GenomeAnalysisTK.jar -T IndelRealigner -R $refgenome -I $inputfolder\/{2}\_{1}_{3}\.bam -targetIntervals $outputfolder\/{2}\_{1}\_intervals\.list -o $outputfolder\/{2}\_{1}\_realigned.bam ::: @samplesnames ::: @pops ::: @term";

}
else{
$cmd2="parallel -j $lnc --results ./logfilesIR --noswap samtools index $inputfolder\/{1}_{2}\.bam ::: @samplesnames ::: @term";
#$cmd3="parallel -j $ncp --results ./logfilesIR --noswap samtools index $inputfolder\/sort_{1}\.bam ::: @samplesnames";
$cmd1 ="parallel -j $ncp --results ./logfilesIR --noswap java -jar ./SNPcallPipe/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $refgenome -I $inputfolder\/{1}_{2}\.bam -o $outputfolder\/{1}\_intervals\.list ::: @samplesnames ::: @term";
$cmd ="parallel -j $ncp --results ./logfilesIR --noswap java -jar ./SNPcallPipe/GenomeAnalysisTK.jar -T IndelRealigner -R $refgenome -I $inputfolder\/{1}_{2}\.bam -targetIntervals $outputfolder\/{1}\_intervals\.list -o $outputfolder\/{1}\_realigned.bam ::: @samplesnames ::: @term";
}
if ($ind=~ m/y/){
print "$cmd2\n";
system ($cmd2);
#system ($cmd3);
}
print "$cmd1\n";
print "$cmd\n";
system ($cmd1);
system ($cmd);
}
1;
