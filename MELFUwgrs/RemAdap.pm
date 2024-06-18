package RemAdap;

my $LEVEL = 1;
sub trim{
my @arg = @_;
foreach $ar (@arg){
        if ($ar =~ /^-i/){our $inf = (split(/ /,$ar))[1];}
        elsif($ar=~ /^-o/){our $outf = (split(/ /, $ar))[1];}
        elsif ($ar=~ /^-fm/){our $fm =(split(/ /, $ar))[1];}
        elsif ($ar=~ /^-nc/){our $nc=(split(/ /, $ar))[1];}
	elsif ($ar=~ /^-exf/){our $exf=(split(/ /, $ar))[1];}
}
#if ($_=~ /^-i$/){$inf=shift @ARGV;}
#elsif ($_=~ /^-o$/){$outf=shift @ARGV;}
#elsif ($_=~ /^-fm$/){$fm=shift @ARGV;}
#elsif ($_=~ /^-nc$/){$nc=shift @ARGV;}
#}
if (not defined ($inf)){print "\nThis script will trimme quality and colapse overlapping PE reads from several samples using AdapterRemoval in parallel. Requires your files after demultiplexing and compressed (sample01.1.fq.gz or sample01_1.fq.gz), all save in one folder.\n\nUsages: RemAdapDTPE.pl\n\t-i <path to inputfolder, if samples are each in one folder use the option -fm>\n\nOptional:\n\t-o <path to outputfolder, default same as inputfolder>\n\t-nc <number the cores or samples to use in parallel, default 4>\n\t-fm <if sequnece file are in one folder per sample (y or n), default n, either the folders or the sequnces shoud be in one folder>\nFor example:\nRemAdapDTPE -i /yuma/WGS/ -o /yuma/remadap/ -fm y -nc 62\n\n";exit;}
print "This is the input folder $inf the output folder $outf the cores $nc and the intruction to move or not $fm\n";
if (not defined ($outf)){$outf=$inf;}
if ( !-d $outf ) {
        `mkdir $outf`;
}
our $tmpdir=$outf."/tmp";
`mkdir $tmpdir`;
if (not defined ($nc)){$nc = 4;}
if (not defined ($fm)){$fm = "n";}
if (not defined ($exf)){$exf = "1.fq.gz,2.fq.gz";}
if ($fm eq "n"){
#my @names=`ls $inf/*1.fq.gz | sed 's/_1.fq.gz//g' | sed 's/.1.fq.gz//g'`;
my @cod=split /\,/, $exf;
#my $test=scalar @cod; 
#print "$cod[0]\t$cod[1]\tYUMA\t$test)\n";
my @names=`ls $inf/*$cod[0] | sed 's/$cod[0]//g'`;
foreach $name (@names){chomp $name; $name=~ s/$inf\///g; push (@nms, $name);}
if ((scalar (@cod)) < 2){
my $cmd="parallel -j $nc --results $tmpdir --tmpdir $tmpdir AdapterRemoval --file1 $inf/{1}$cod[0] --basename $outf/{1} --minlength 70 --trimns --trimqualities --collaps --gzip ::: @nms";
print "$cmd\n\n";
`parallel -j $nc --results $tmpdir --tmpdir $tmpdir AdapterRemoval --file1 $inf/{1}$cod[0] --basename $outf/{1} --minlength 70 --trimns --trimqualities --collapse --gzip ::: @nms`;
### yuma
}
elsif ((scalar (@cod)) == 2){
my $cmd="parallel -j $nc --results $tmpdir --tmp $tmpdir AdapterRemoval --file1 $inf/{1}$cod[0] --file2 $inf/{1}$cod[1] --basename $outf/{1} --minlength 70 --trimns --trimqualities --collapse --gzip ::: @nms";
print "$cmd\n\n";
`parallel -j $nc --results $tmpdir --tmpdir $tmpdir AdapterRemoval --file1 $inf/{1}$cod[0] --file2 $inf/{1}$cod[1] --basename $outf/{1} --minlength 70 --trimns --trimqualities --collapse --gzip ::: @nms`;
}
}
else {
my @cod=split /\,/, $exf;
my @names=`ls -d $inf/*/`;

foreach $name (@names){chomp $name; $name=~ s/$inf\///g; $name=~ s/\///g; push (@nms, $name);  
	my $filesN = `ls $inf/$name/* | wc -l`;
        if ($filesN <= 3){`mv $inf/$name/$name\_*_1.fq.gz $inf/$name/$name\_1.fq.gz`; `mv $inf/$name/$name\_*_2.fq.gz $inf/$name/$name\_2.fq.gz`;}
	else {`cat $inf\/$name\/$name\_\*_1.fq.gz \>\> $inf\/$name\/$name\_1.fq.gz`; `cat $inf\/$name\/$name\_\*_2.fq.gz \>\> $inf\/$name\/$name\_2.fq.gz`;}
}

my $cmd="parallel -j $nc --results $tmpdir --tmp $tmpdir AdapterRemoval --file1 $inf/{1}/{1}$cod[0] --file2 $inf/{1}$cod[1] --basename $outf/{1} --minlength 70 --trimns --trimqualities --collapse --gzip ::: @nms";
print "$cmd\n\n";
`parallel -j $nc --results $tmpdir --tmpdir $tmpdir AdapterRemoval --file1 $inf/{1}/{1}$cod[0] --file2 $inf/{1}/{1}$cod[1] --basename $outf/{1} --minlength 70 --trimns --trimqualities --collapse --gzip ::: @nms`;
}


} 

1;
