package bedmarkrep;

my $LEVEL = 1;
sub markrep{
my @arg = @_;
foreach $ar (@arg){
	if ($ar =~ /^-i /){our $inf = (split(/ /,$ar))[1];}
	elsif($ar=~ /^-o /){our $outf = (split(/ /, $ar))[1];}
	elsif($ar=~ /^-rg /){our $refgr = (split(/ /, $ar))[1];}
	elsif ($ar=~ /^-snc /){$snc=(split(/ /, $ar))[1];}
}

if (not defined ($inf && $refgr)){print "\nThis script will markrepeat regions using bedtools in parallel. Requires a repeat masked reference genome in a folder; and the samples  bam files after mapping, merge and mark duplicates all save in one folder.\n\nUsages: bedtoolsrepeats.pl\n\t-i <path to inputfolder>\n\t-rg <path to reference repeats masked>\n\nOptional:\n\t-o <path to outputfolder, default bedout>\n\t-nc <number the cores or samples to use in parallel, default 4>\n\nFor example:\n\nbedtoolsrepeats.pl -i /yuma/dedup/ -o /yuma/repeatmrk/ -rg /yuma/makosharkgenome/repeat_mask_Iocy.bed -snc 60\n\n";exit;}
if (not defined ($outf)){$outf="./bedout";}
our $tmpdir=$outf."/tmp";
#print "$tmpdir\n";
if ( !-d $outf ) {
	`mkdir $outf`;
}
if ( !-d $tmpdir ) {
	`mkdir $tmpdir`;
}
if (not defined ($snc)){$snc=4;}
my @names=`ls $inf/*.bam | sed 's/__/_/g' | sed 's/_markdup.bam//g' | sed 's/_mergeA.bam//g'`;
foreach $name (@names){chomp $name; $name=~ s/$inf\///g; push (@nms, $name);}
my $cmd="parallel -j $nc bedtools intersect -abam $inf/{1}*.bam -b $refgr -v \'\>\' $outf/{1}_repeatmasked.bam ::: @nms";
`parallel -j $snc bedtools intersect -abam $inf/{1}*.bam -b $refgr -v \'\>\' $outf/{1}_repeatmasked.bam ::: @nms`;
`parallel -j $snc samtools index $outf/{1}_repeatmasked.bam.bam ::: @nms`;
}
1;
