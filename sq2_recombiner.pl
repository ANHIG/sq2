#!/usr/bin/perl
#
# A perl script to identify recombinant regions as described in;
#
# http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006862
#
# Please see README and LICENSE for details of use
#
#=========================================================================================

# Example input file is a_recombinant.csv which is a csv file of allele name and postion,
# the first line is a line details the postions included. 

use List::Util 'max';
use List::Util 'min';
use List::Util 'sum';
use List::Util 'first';
use Array::Compare;
use Array::Utils qw(:all);
use Array::Unique;
use Data::Dumper;

use sort_new_nomenc;

$lookup = @ARGV[0];
chomp($lookup);

$junk_data = 1;
open(IN, $lookup);
while(<IN>)
  {  
  #chomp();
  if ( $_ =~ /Recombinant profile/ )
    {
    print "$_\n";
    if ( $tracts =~ /[123]/ && $junk_data eq 0 )
      {
      $recomb_forward{$allele} = $tracts;
      }
    $tracts = '';
    $junk_data = 0;
    $allele = substr($_,24,);
    $allele =~ s/ .*//;
    chomp($allele);
    }
  elsif ( $_ =~ /\t[123]\t[ACTG].*\t\d+\.\.\d+/ )
    {
    #print "$_\n";
    #print "GU = @genes_used\n";
    $_ =~ s/\tHLA.*//;
    $tracts .= $_;
    }
  elsif ( $_ =~ /\t4\t[ACTG].*\t\d+\.\.\d+/ || $_ =~ /unique/ )
    {
    $tracts = '';  
    $junk_data = 1;
    }
  }
if ( $tracts =~ /[123]/ )
  {
  $recomb_forward{$allele} = $tracts;
  }  

$lookup =~ s/forward/reverse/;
$junk_data = 1;
open(IN, $lookup);
while(<IN>)
  {  
  if ( $_ =~ /Recombinant profile/ )
    {
    if ( $tracts =~ /[123]/ )
      {
      $recomb_reverse{$allele} = $tracts;
      }
    $tracts = '';
    $junk_data = 0;
    $allele = substr($_,24,);
    $allele =~ s/ .*//;
    chomp($allele);
    }
  elsif ( $_ =~ /\t[123]\t[ACTG].*\t\d+\.\.\d+/ )
    {
    $_ =~ s/\tHLA.*//;
    $tracts .= $_;
    }
  elsif ( $_ =~ /\t4\t[ACTG].*\t\d+\.\.\d+/ || $_ =~ /unique/ )
    {
    $tracts = '';  
    $junk_data = 1;
    }    
  }
if ( $tracts =~ /[123]/ )
  {
  $recomb_reverse{$allele} = $tracts;
  }  
  
$version_0_count = 0;  
$version_12_count = 0;  
$version_123_count = 0;
  
foreach $allele ( %recomb_forward )
  {
  if ( $allele =~ /HLA/ )
    {    
    if ( $recomb_forward{$allele} =~ /3\t[ACTG]/ )
      {
      ++$version_123_count;
       #print "Version 123 allele\n";
      print "Version 123 allele\n";
      print "$allele\n";
      print "FORWARD\t$recomb_forward{$allele}\n";
      print "REVERSE\t$recomb_reverse{$allele}\n";
      push(@version123,$allele);

      #print "\Min Bounds is:";

      $min_start = $recomb_forward{$allele};
      $min_start =~ s/[\n\r]//g;
      $min_start =~ s/.*\t2\t[ACTG]{1,}\t//;
      $min_start =~ s/\.\..*//;
      
      $min_end = $recomb_reverse{$allele};
      $min_end =~ s/[\n\r]//g;
      $min_end =~ s/.*\t2\t[ACTG]{1,}\t//;
      $min_end =~ s/\.\..*//;
      
      $min_size = $min_end-$min_start;
      ++$min_size;
      
      if ( $min_start eq $min_end || $min_size < 1 )
        {
        print "$allele ($min_start eq $min_end )\n";
        #exit;
        }
      
      print "Min: $min_start to $min_end ($min_size)\n";
      for ( $ms=76; $ms<619; ++$ms )
        {
        if ( $ms > $min_start && $ms < $min_end )
          {
          $gene_conv_min{$ms}++;
          }
        }
      push(@min_size, $min_size);
      
      $max_end = substr($recomb_forward{$allele}, -9, );
      $max_end =~ s/[ACGT\t\n\r]//g;
      $max_end =~ s/\.\..*//;
      --$max_end;
      $max_start = substr($recomb_reverse{$allele}, -9, );
      $max_start =~ s/[ACGT\t\n\r]//g;
      $max_start =~ s/\.\..*//;  
      ++$max_start;
      
      $max_size = $max_end-$max_start; 
      ++$max_size;
      print "Max: $max_start to $max_end ($max_size)\n";
       for ( $ms=76; $ms<619; ++$ms )
        {
        if ( $ms > $max_start && $ms < $max_end )
          {
          $gene_conv_max{$ms}++;
          }
        }
      push(@max_size, $max_size);
 	  #exit;
      print "XLS2\t$allele\t76\t$max_start\t$min_start\t$min_end\t$max_end\t618\n";
      
      print "XLS\t$allele\t$max_start/$min_start\t$min_end/$max_end\n";
      }
    elsif ( $recomb_forward{$allele} =~ /2\t[ACTG]/ )
      {
      ++$version_12_count;
      print "Version 12 allele\n";
      print "$allele\n";
      print "FORWARD\t$recomb_forward{$allele}\n";
      print "REVERSE\t$recomb_reverse{$allele}\n";
        
      push(@version12,$allele);

      #print "\tMax Bounds is:";

      $max_end = substr($recomb_forward{$allele}, -9, );
      $max_end =~ s/[ACGT\t\n\r]//g;
      $max_end =~ s/\.\..*//;
      $max_start = substr($recomb_reverse{$allele}, -9, );
      $max_start =~ s/[ACGT\t\n\r]//g;
      $max_start =~ s/\.\..*//;
            
      for ( $ms=76; $ms<619; ++$ms )
        {
        if ( $ms > $max_start && $ms < $max_end )
          {
          $single_recomb_hits{$ms}++;
          }
        }
      }  
     else
      {
      ++$version_0_count;
      }
    }
  }  
  
print "Version 0 count is $version_0_count\n";  
print "Version 12 count is $version_12_count\n";  
print "Version 123 count is $version_123_count\n";  

print "VERSION12\t@version12\n";
print "VERSION123\t@version123\n";
exit;

# last;
# 
 print "Data for single recombinant histograms\n";
 
 for ( $ms=76; $ms<619; ++$ms )
   {
   print "$ms\t $single_recomb_hits{$ms}\n";  
   }
 
 print "Data for minimum double recombinant histograms\n";#
 
 for ( $ms=76; $ms<619; ++$ms )
   {
   print "$ms\t$gene_conv_min{$ms}\n";  
   }
   
 print "Data for maximum double recombinant histograms\n";#
 
 for ( $ms=76; $ms<619; ++$ms )
   {
   print "$ms\t$gene_conv_max{$ms}\n";  
   }

 print "Data for Figure 11 histograms\n";#
 
 for ( $ms=76; $ms<619; ++$ms )
   {
   print "$ms\t$gene_conv_max{$ms}\t$gene_conv_min{$ms}\n";  
   }

  
$min_min = min @min_size;
$max_min = max @min_size;
$avg_min = sum(@min_size)/@min_size;

print "Min Size: @min_size\n";
print "Min Size: Min $min_min; Max $max_min; Avg $avg_min\n";

$min_max = min @max_size;
$max_max = max @max_size;
$avg_max = sum(@max_size)/@max_size;

print "Max Size: @max_size\n";
print "Max Size: Min $min_max; Max $max_max; Avg $avg_max\n";

      