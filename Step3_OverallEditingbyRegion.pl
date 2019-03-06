#!/usr/bin/perl -w
#use strict;
my($six,$twelve,%hash_up_v2,%hash_down_v2,$ID,$total,$out, $count, $sum, $filecount,$outread,$n,$p1,$sqtotal,$editreads,$editlevel,$annot, $std,$p2, $p3, $p4, $p5, $p6,$p7,$p8,$p9,$p10,$poor,$cov, $length, $GC,%hash,%counter,$control, $pval,$min, $name,%hash_up, %hash_down, @Pvals, @sorted_pvals, $line, $gene, $beta, $t,$p,$gi5,$file);

my $count3UTR=0;
my $sum3UTR =0;
my $sum_editreads3UTR=0;

my $count5UTR=0;
my $sum5UTR =0;
my $sum_editreads5UTR=0;

my $countExon=0;
my $sumExon =0;
my $sum_editreadsExon=0;

my $countIntron=0;
my $sumIntron =0;
my $sum_editreadsIntron=0;

my $countIntergenic=0;
my $sumIntergenic =0;
my $sum_editreadsIntergenic=0;

my $countncRNA=0;
my $sumncRNA =0;
my $sum_editreadsncRNA=0;

my $overallediting3UTR =0;
my $overallediting5UTR =0;
my $overalleditingExon =0;
my $overalleditingIntron =0;
my $overalleditingIntergenic =0;
my $overalleditingncRNA =0;

open(INFILE, $ARGV[0]) or die"File1 is Dead\n";
$file=$ARGV[0];
while(<INFILE>)
{
            $line=$_;
            chomp $line;
            if ($line !~m/(N\/A)$/)
            {
            ($gi)= ($line=~m/\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/); #print "\n$gi\n"; die;
            ($cov)=($line=~m/(\S+)\s+\S+\s+\S+$/);
            ($editreads)=($line=~m/(\S+)\s+\S+$/);
            ($editlevel)=($line=~m/(\S+)$/);
            ($annot)=($line=~m/\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);
                        
                        if ($gi=~m/3UTR/ and $line !~m/(N\/A)$/ and $line !~m/#chrom/)
                        {
                        $sum3UTR+=$cov; #print "\n$line\n"; die;
                        $sum_editreads3UTR+=$editreads;
                        $count3UTR++;
                        }
                        
                        if ($gi=~m/5UTR/ and $line !~m/(N\/A)$/ and $line !~m/#chrom/)
                        {
                        $sum5UTR+=$cov;
                        $sum_editreads5UTR+=$editreads;
                        $count5UTR++;
                        }

                        if ($gi=~m/Syn/ or $gi=~m/Nonsyn/ or $gi=~m/exonic/ and $line !~m/(N\/A)$/ and $line !~m/#chrom/)
                        {
                        $sumExon+=$cov;
                        $sum_editreadsExon+=$editreads;
                        $countExon++; }
                        }
            
                        if ($gi=~m/intronic/ and $line !~m/(N\/A)$/ and $line !~m/#chrom/)
                        {
                        $sumIntronic+=$cov;
                        $sum_editreadsIntronic+=$editreads;
                        $countIntronic++; 
                        }

                        if ($gi=~m/intergenic/  and $line !~m/(N\/A)$/ and $line !~m/#chrom/)
                        {
                        $sumIntergenic+=$cov;
                        $sum_editreadsIntergenic+=$editreads;
                        $countIntergenic++; 
                        }

                        if ($gi=~m/ncRNA/  and $line !~m/(N\/A)$/ and $line !~m/#chrom/)
                        {
                        $sumncRNA+=$cov;
                        $sum_editreadsncRNA+=$editreads;
                        $countncRNA++; 
                        }
            
}

#If coverage is poor, print the sample name... could be an outlier sample and requires individual inspection
if ($count3UTR <10)     {print  "$file\tPoor coverage\n";}
if ($count5UTR <10)     {print  "$file\tPoor coverage\n";}
if ($countExon <10)     {print  "$file\tPoor coverage\n";}
if ($countIntronic <10) {print  "$file\tPoor coverage\n";}
if ($countIntergenic <10){print  "$file\tPoor coverage\n";}
if ($countncRNA <10)    {print  "$file\tPoor coverage\n";}


else{
     eval{  #compute overall editing rates
            $overallediting3UTR = $sum_editreads3UTR/$sum3UTR;
            $overallediting5UTR = $sum_editreads5UTR/$sum5UTR;
            $overalleditingExon = $sum_editreadsExon/$sumExon;
            $overalleditingIntron = $sum_editreadsIntronic/$sumIntronic;
            $overalleditingIntergenic = $sum_editreadsIntergenic/$sumIntergenic;
            $overalleditingncRNA = $sum_editreadsncRNA/$sumncRNA;
 
#print out file name, total detected sites, total number of edited reads, total coverage, overall editing rates
print  "$file\t$count3UTR\t$sum_editreads3UTR\t$sum3UTR\t$overallediting3UTR\t3UTR\n";
print  "$file\t$count5UTR\t$sum_editreads5UTR\t$sum5UTR\t$overallediting5UTR\t5UTR\n";
print  "$file\t$countExon\t$sum_editreadsExon\t$sumExon\t$overalleditingExon\tExonic\n";
print  "$file\t$countIntronic\t$sum_editreadsIntronic\t$sumIntronic\t$overalleditingIntron\tIntronic\n";
print  "$file\t$countIntergenic\t$sum_editreadsIntergenic\t$sumIntergenic\t$overalleditingIntergenic\tIntergenic\n";
print  "$file\t$countncRNA\t$sum_editreadsncRNA\t$sumncRNA\t$overalleditingncRNA\tncRNA\n";

}
}





