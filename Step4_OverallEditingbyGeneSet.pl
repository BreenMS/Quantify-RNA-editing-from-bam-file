#!/usr/bin/perl -w
#use strict;
my($six,$twelve,%hash_up_v2,%hash_down_v2,$ID,$total,$out, $count, $sum, $filecount,$outread,$n,$p1,$sqtotal,$editreads,$editlevel,$annot, $std,$p2, $p3, $p4, $p5, $p6,$p7,$p8,$p9,$p10,$poor,$cov, $length, $GC,%hash,%counter,$control, $pval,$min, $name,%hash_up, %hash_down, @Pvals, @sorted_pvals, $line, $gene, $beta, $t,$p,$gi5,$file);

open(INFILE, $ARGV[0]) or die"File1 is Dead\n";
while(<INFILE>)
{
$line=$_;
chomp $line;
if ($line=~m/\S+/)
{
    ($gi)=($line=~ m/(\S+)/);
    ($symbol)=($line=~ m/\S+\s+(\S+)/);
    $hash{$gi}=$symbol;
}
}

my $countG=0;
my $sumG =0;
my $sum_editreadsG=0;
my $sum_editlevelsG=0;

my $countS=0;
my $sumS =0;
my $sum_editreadsS=0;
my $sum_editlevelsS=0;

my $countD=0;
my $sumD =0;
my $sum_editreadsD=0;
my $sum_editlevelsD=0;

my $overalleditingG=0;
my $overalleditingD=0;
my $overalleditingS=0;

open(INFILE1, $ARGV[1]) or die"File1 is Dead\n";
$file=$ARGV[1];
while(<INFILE1>)
{
            $line=$_;
            chomp $line;
            if ($line !~m/(N\/A)$/)
            {
            ($gi)= ($line=~m/\S+\s+\S+\s+(\S+)/); 
            ($cov)=($line=~m/(\S+)\s+\S+\s+\S+$/);

            ($editreads)=($line=~m/(\S+)\s+\S+$/);
            ($editlevel)=($line=~m/(\S+)$/);
            ($annot)=($line=~m/\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/);

                        if ($hash{$gi} and $line !~m/(N\/A)$/ and $line !~m/#chrom/)
                        {
                        if($hash{$gi} eq "Glutamatergic"){
                        $sumG+=$cov;
                        $sum_editreadsG+=$editreads;
                        #$sum_editlevelsG+=$editlevel;
                        $countG++; }
                        
                        if($hash{$gi} eq "Serotonergic"){
                        $sumS+=$cov;
                        $sum_editreadsS+=$editreads;
                        #$sum_editlevelsS+=$editlevel;
                        $countS++; }

                        if($hash{$gi} eq "Dopaminergic"){
                        $sumD+=$cov;
                        $sum_editreadsD+=$editreads;
                        #$sum_editlevelsD+=$editlevel;
                        $countD++; }
                        }
            }
}

#If coverage is poor, print the sample name... could be an outlier sample and requires individual inspection
if ($countG <10){print  "$file\tPoor coverage\n";}
if ($countS <10){print  "$file\tPoor coverage\n";}
if ($countD <10){print  "$file\tPoor coverage\n";}

else{
     eval{  #compute overall editing rates
            $overalleditingG = $sum_editreadsG/$sumG;
            $overalleditingS = $sum_editreadsS/$sumS;
            $overalleditingD = $sum_editreadsD/$sumD;

#print out file name, total detected sites, total number of edited reads, total reads, overall editing rates (total edited reads/total reads)
          print  "$file\t$countG\t$sum_editreadsG\t$sumG\t$overalleditingG\tGlutamatergic\n";
          print  "$file\t$countS\t$sum_editreadsS\t$sumS\t$overalleditingS\tSerotonergic\n";
          print  "$file\t$countD\t$sum_editreadsD\t$sumD\t$overalleditingD\tDopaminergic\n";}
    }

