#!/usr/bin/perl -w
#use strict;

my($line, $file, $total_sites,$sum_coverage,$sum_editreads, $sum_editlevels,$overallediting);

$total_sites=0;
$sum_coverage =0;
$sum_editreads=0;
$sum_editlevels=0;
$overallediting=0;

print "SampleName\tTotalDetectedSites\tTotalEditedReads\tTotalCoverage\tOverallEditing\tCategory\n";   

open(INFILE, $ARGV[0]) or die"File1 is Dead\n";
$file=$ARGV[0];
while(<INFILE>)
{
$line=$_;
chomp $line;

            ($coverage)=($line=~m/(\S+)\s+\S+\s+\S+$/);
            ($editreads)=($line=~m/(\S+)\s+\S+$/);

                        if ($line !~m/(N\/A)$/ and $line !~m/#chrom/)
                        {
                        $sum_coverage+=$coverage;
                        $sum_editreads+=$editreads;
                        #$sum_editlevels+=$editlevel;
                        $total_sites++; 
                        }          

}

if ($total_sites <10) 
{print "$file\tPoor coverage...likely a sample outlier\n";}

else{
            #compute overall editing rates
            $overallediting = $sum_editreads/$sum_coverage; 
            
            #print out file name, total detected sites, total number of edited reads, total coverage, overall editing rates
            print "$file\t$total_sites\t$sum_editreads\t$sum_coverage\t$overallediting\tOverallEditingRates\n";   
    }
