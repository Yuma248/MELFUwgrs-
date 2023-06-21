#!~/miniconda3/envs/SNPcallPipeY/bin/perl -w
while (@ARGV){
	$_=shift @ARGV;
	if ($_=~ /^-stp$/){$stp=shift @ARGV;}
	elsif ($_=~ /^-i$/){$input=shift @ARGV;}
	elsif($_=~ /^-o$/){$output=shift @ARGV;}
	elsif ($_=~ /^-fm$/){$fm=shift @ARGV;}
	elsif ($_=~ /^-lnc$/){$lnc=shift @ARGV;}
	elsif ($_=~ /^-snc$/){$snc=shift @ARGV;}
	elsif ($_=~ /^-rg$/){$reference=shift @ARGV;}
	elsif ($_=~ /^-stpn$/){$stpn=shift @ARGV;}
	elsif ($_=~ /^-pf$/){$pf=shift @ARGV;}
	elsif ($_=~ /^-al$/){$B=shift @ARGV;}
	elsif ($_=~ /^-t$/){$type=shift @ARGV;}
	elsif ($_=~ /^-ind$/){$ind=shift @ARGV;}
	elsif ($_=~ /^-exf$/){$exf=shift @ARGV;}
	elsif ($_=~ /^-chrf$/){$chrf=shift @ARGV;}
}
if (not defined ($stp)){print "\nThis pipeline  trim and\/or map reads, and call SNPs/SVs and filter them. The required arguments and inputs depend of the steps you want to perform. You will have to use the option -stpn and select the step you want to run. To check the step names run the script without arguments, if you want to check the arguments for one step just run the script with the specific step but without extra arguments (example MELFUwgrs.pl -stp trim).\n\nUsage:\nMELFUwgrs.pl\n\t-stp <You need at least determine what steps you want to run>\n\t\tindref\: <Indexs the reference genome with samtools, picard, bowtie2 and snap>\n\t\ttrim\: <It will use AdapterRemoval to trim and filter reads>\n\t\tconcat\: <It will concatenate fastq files of the same sample but different runs in one file>\n\t\talignment\: <It will use bowtie2, bwa or snap to align reads to a reference genome>\n\t\tdedup\: <This step will sort sam/bam files, convert sam to bam (if necessary) and mask duplicates>\n\t\tindelrea\: <This step will locally realign indels, although this is not recommended any more>\n\t\tbedmarkrep: <This step will mask repeat regions in the genome>\n\t\tsnpcalling\: <This step will use ANGSD to simultaneously call and genotype SNP>\n\t\tfiltering\: <This step will use vcftools to filter SNPs, I recommend to use this automatically to have an idea of your data, but play whit the parameters if you have the time>\n\n"; exit;}
if (not defined ($stprn)){$stprn = 0};
if (not defined ($stpn)){$stpn =0};
use File::Basename;
use lib dirname (__FILE__) . "/MELFUwgrs";
use Cwd 'abs_path';
my $SCP1= abs_path($0);
$SCP1 =~ s/\.pl$/\//;
our $stprn =0;
our @stptr = split (/,/, $stp);
foreach $stp (@stptr){
	if ($stp eq "indref"){
        	use indexgenome;
		if (not defined ($input && $reference)){print "\nThis script will index a reference genome using samtools, picard, bowtie2 and snap. The genome should have extension fna.\n\nUsage:\nMELFUwgrs.pl -stp indref\n\t-i <input folder where the reference is, and where all the output indexes, dictionaries and folders will be saved>\n\t-rg <the name of the reference genome including the extension>\n\t-pf <path to the picar jar executable if different to MELFUwgrs>\n\nExample:\nMELFUwgrs.pl -stp indref -i yuma/genomes/ -rg Taust.Chr.fna -pf /local/SNOcallPipe/\n"; exit;}
		if (not defined ($pf)){$pf = $SCP1;}
        	our @arg = ("-i $input","-rg $reference","-pf $pf");
	        indexgenome::indrg(@arg);
        	$stprn = 0;
	}
	elsif ($stp eq "trim" or $stprn == 1){
		use RemAdap;
		if (not defined ($input)){print "\nThis script will trim quality and collapse overlapping PE reads from several samples using AdapterRemoval in parallel. Requires your files compressed (sample01.1.fq.gz or sample01_1.fq.gz), all save in one folder, or each sample in a folders and all the folder per sample in one parental folder.\n\nUsages:\nMELFUwgrs.pl -stp trim\n\t-i <path to inputfolder, if samples are each in one folder use the option -fm>\n\nOptional:\n\t-o <path to outputfolder, default same as inputfolder>\n\t-snc <number the cores or samples to use in parallel, default 4>\n\t-fm <if sequence files are in one folder per sample (y or n), default n, either the folders or the sequences should be in one parental folder>\n\t-exf <Extension format of sequences fastq files, default 1.fq.gz,2.fq.gz>\n\nExample:\nMELFUwgrs.pl -stp trim -i /yuma/WGS/ -o /yuma/remadap/ -fm y -snc 62 -exf F.fq.gz,R.fq.gz\n\n";exit;}
		if (not defined ($output)){$output=$input;}
		if (not defined ($snc)){$snc = 4;}
		if (not defined ($fm)){$fm = "n";}
		if (not defined ($exf)){$exf ="1.fq.gz,2.fq.gz";}
        	our @arg = ("-i $input","-o $output","-fm $fm","-nc $snc","-exf $exf");
	        RemAdap::trim(@arg);
        	$stprn = 2;
	}
	elsif ($stp eq "concat" or $stprn == 1){
		use concat;
		if (not defined ($inputfolder && $outputfolder)){print "\nThis script will concatenate fastq files of raw or trimmed reads of several samples in parallel. It requires fastq files of all samples stored in the same folder\n\nUsage:\nMELFUwgrs.pl -stp concat\n\t-i <path to the folder with the input fastq files>\n\t-o <path to the output folder>\n Optional:\n\t-snc <number of runs in parallel, default 10>\n\t-t <method used for trimming. Trimmomatic TR, AdapterRemoval AR or None NO if are raw sequences, default AR>\n\t-exf <this will tell the script how to extract the name of each sample, and should include all extra information at the end of the file names that is not related to the sample name, default P1_L001_>\n\nFor example:\nMELFUwgrs.pl -stp concat -i /home/Yumafan/demultiplex/trimmed/ -o /home/Yumafan/concatenated -snc 10 -t AR -exf P1_L001_\n\n"; exit;}
		if (not defined ($snc)){$snc =10;}
		if (not defined ($type)){$type = "RA";}
		if (not defined ($exf)){$exf ="P1_L001_";}
        	our @arg = ("-i $input","-o $output","-nc $snc","-exf $exf","-t $type");
	       concat::coca(@arg);
        	$stprn = 2;
	}
	elsif ($stp eq "mergbam" or $stprn == 1){
                use mergbam;
		if (not defined ($inputfolder && $outputfolder)){print "\nThis script will merge bam files of different lanes or runs of the same sample, it is important that you check the runs have similar quality and have the same category of read (see samtools flagstat). The script can work with several samples in parallel. It requires bam files of all samples stored in the same folder\n\nUsage:\nMELFUwgrs.pl -stp merge\n\t-i <path to the folder with the input bam files>\n\t-o <path to the output folder>\n Optional:\n\t-snc <number of runs in parallel, default 10>\n\t-exf <this will tell the script how to extract the name of each sample, and should include all extra information at the end of the file names that is not related to the sample name, default _S*P1_L001_>\n\t-ind <if you need to samtools index the input files, default yes>\n\nFor example:\nMELFUwgrs.pl -stp concat  -i /home/Yumafan/dedupout/ -o /home/Yumafan/mergedbam -snc 10 -exf S*P1_L001__markdup -ind yes\n\n"; exit;}
                if (not defined ($snc)){$snc =10;}
                if (not defined ($exf)){$exf ="_S\*P1_L001__markdup";}
                if (not defined ($ind)){$ind ="yes";}
                our @arg = ("-i $input","-o $output","-snc $snc","-exf $exf", "-ind $ind");
                mergbam::merg(@arg);
                $stprn = 2;
	}
	elsif ($stp eq "alignment" or $stprn == 2){
        	use BWAB2;
		if (not defined ($reference && $input && $output)){print "\nThis script will map reads of several samples to a reference genome using bwa-mem, bowtie2 or snap-aligner in parallel. It requires files (after trimming) of all samples stored in the same folder\n\nUsage:\nMELFUwgrs.pl -stp alignment\n\t-rg <path to the reference genome fasta file, or reference genome folder if SNAP is used>\n\t-i <path to the folder with the input fasta files>\n\t-o <path to the output folder>\n Optional:\n\t-lnc <number of runs in parallel, default 1>\n\t-snc <number cores for each run, default 4. The total number of cpus to be used will be lnc * snc>\n\t-al <aligner to be used, BWA,  B2 (bowtie2), or SNAP, default BWA>\n\t-t <sequencing type single-end \"S\" or paired-end \"P\", default P> \n\nFor example:\nMELFUwgrs.pl -stp alignment -rg /home/Yumafan/genome/reference_genome.fasta -i /home/Yumafan/demultiplex/trimmed/ -o /home/Yumafan/bwaout/ -lnc 10 -snc 4 -al BWA -t P\n\n"; exit;}
		if (not defined ($lnc)){$lnc=1;}
		if (not defined ($snc)){$snc=4;}
		if (not defined ($B)){$B="BWA";}
		if (not defined ($type)){$type="P";}
	        our @arg = ("-i $input","-o $output","-rg $reference","-al $B","-t $type","-ncr $snc","-ncp $lnc");
        	BWAB2::align(@arg);
	        $stprn = 3;

	}
	elsif ($stp eq "dedup" or $stprn == 3){
        	use samdedup;
		if (not defined ($input)){print "\nThis script will convert from sam to bam (if required), sort by name, fix mates, sort by coordinates and mark duplicates using samtools in parallel. Requires your sam/bam files after mapping (sample01.sam/bam), all save in one folder.\n\nUsages:\nMELFUwgrs.pl -stp dedup\n\t-i <path to inputfolder>\n\nOptional:\n\t-o <path to outputfolder, default samout>\n\t-lnc <number the cores or samples to use in parallel, default 4>\n\t-snc <number of cores per sample, default 1. The total number of cpus to be used will be lnc * snc >\n\nExample:\nMELFUwgrs.pl -stp dedup -i /home/Yumafan/demultiplex/trimmed/bwaout -o /home/Yumafan/dedupout/ -lnc 10 -snc 4\n\n";exit;}
		if (not defined ($lnc)){$lnc=4;}
		if (not defined ($snc)){$snc=1;}
		if (not defined ($output)){$output="./samout";}
	        our @arg = ("-i $input","-o $output","-lnc $lnc","-snc $snc");
        	samdedup::dedup(@arg);
	        $stprn = 4;

	}
	elsif ($stp eq "indelrea" or $stprn == 4){
        	use indelrealigment;
		if (not defined ($reference && $input && $output)){print "\nThis script will create a list of indels per sample, and it will realign reads around indels for a more efficient variants calling. It uses GATK in parallel to process several samples at the same time, this step is not recommend anymore.\n\nUsage:\nMELFUwgrs.pl -stp indelrea\n\t-i <input folder with sorted bam files>\n\t-o <output folder to save the realignments>\n\t-rg <path to the reference genome>\nOptional:\n\t-snc <number cores to run in parallel, default 4>\n\t-lnc <Some of the processes are heavy in RAM usage which limits the number of runs in parallel. But you can speed the first part using more cores with this parameter, default equal to -snc value. The total number of cpus to be used will be the biggest value of -lnc and -snc>\n\t-ind <\"yes\" if you need to index your bam files, or \"no\" if they are already indexed, default \"yes\">\n\nExample:\nMELFUwgrs.pl -stp indelrea -i ./Yuma/sam2bamout/ -o ./Yuma/indelout/ -rg ./Yuma/genomes/reference.fasta -snc 8 -ind yes\n\n"; exit;}
		if (not defined ($snc)){$snc=8;}
		if (not defined ($lnc)){$lnc=$snc;}
		if (not defined ($ind)){$ind="yes";}	
        	our @arg = ("-i $input","-o $output","-ind $ind","-ncp $snc","-rg $reference", "-lncp $lnc");
	        indelrealigment::indelreal(@arg);
        	$stprn = 5;
	}
        elsif ($stp eq "bedmarkrep" or $stprn == 5){
                use bedmarkrep;
		if (not defined ($input && $reference)){print "\nThis script will mark repeat regions using bedtools in parallel. Requires a repeat masked reference genome in a folder; and the samples  bam files after mapping, merge and mark duplicates all save in one folder.\n\nUsages: bedtoolsrepeats.pl\n\t-i <path to inputfolder>\n\t-rg <path to reference repeat regions masked, BAM/BED/GFF/VCF>\n\nOptional:\n\t-o <path to outputfolder, default bedout>\n\t-nc <number the cores or samples to use in parallel, default 4>\n\nFor example:\n\nbedtoolsrepeats.pl -i /yuma/dedup/ -o /yuma/repeatmrk/ -rg /yuma/makosharkgenome/repeat_mask_Iocy.bed -snc 60\n\n";exit;}
                if (not defined ($snc)){$snc =4;}
                if (not defined ($output)){$output ="./bedout";}
                our @arg = ("-i $input","-o $output","-snc $snc","-rg $reference");
                bedmarkrep::markrep(@arg);
                $stprn = 6;
        }
 	elsif ($stp eq "snpcalling" or $stprn == 6){
        	use SNPangsd;
		if (not defined ($input && $output && $reference)){print "\nYUMA\n\n\n\nThis script will create the ANGSD command to call and genotype SNPs for several samples in parallel, it use a likelihood approach recommended for low coverage genome sequencing. Its need the mapped bam files in a folder,and a indexed reference genome. It is recommended to use  just chromosomes or the biggest scaffolds, so you can use a file with the names of chromosomes or scaffolds to be used, default will use the first 20 chromosome, scaffolds or contigs in the reference genome.\n\nUsage:\nMELFUwgrs -stp calling\n\t-i <input folder with mapped bam files with realigned indels>\n\t-o <output folder to save vcf files>\n\t-rg <reference genome>\n\t-chf <File with chromosome names list, default the first 20 scaffold or contigs in the reference will be used>\n\t-snc <number the cores to be used in parallel, recommend to use the number of Chromosomes, default 20>\n\nExample:\nMELFUwgrs.pl -stp snpcalling -i ./Yuma/indelrealigned/ -o ./Yuma/rawsnp/ -rg ./Yuma/genome/reference_genome.fasta -snc 23 -chrf\n\n"; exit;}
		if (not defined ($snc)){$snc = 20;}
		if (not defined ($chrf)){$chrf = "N";}
	        our @arg = ("-i $input","-o $output","-nc $snc","-rg $reference", "-chrf $chrf");
        	SNPangsd::callSNPs(@arg);
	        $stprn = 7;
	}
	elsif ($stp eq "filtering" or $stprn == 7){
        	##use vcftoolsF;
	        ##our @arg = ("-i $input","-o $output","-ind $ind","-ncp $lnc","-rg $reference");
        	##vcftoolsF::filteringY(@arg);
	}
}
