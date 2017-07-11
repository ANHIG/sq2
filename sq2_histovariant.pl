#!/usr/bin/perl
#
# A perl script to generate variation histograms as described in;
#
# http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006862
#
# Please see README and LICENSE for details of use
#
#=========================================================================================

use List::Util 'max';
use Array::Compare;
use Array::Utils qw(:all);
use Array::Unique;

#use sort_new_nomenc;

$omit         = 0;
$pos_offset   = 74;

foreach $argv (@ARGV) 
  {
  if ( $argv =~ /^i/ )
    {
    $argv =~ s/i=//;    
    $filein = $argv;
    }
  elsif ( $argv =~ /^s/ )
    {
    $argv =~ s/s=//;    
    $remove_snps = $argv;
    }
  elsif ( $argv =~ /^o/ )
    {
    $argv =~ s/o=//;    
    $omit = $argv;
    }  
  elsif ( $argv =~ /^f2/ )
    {
    $argv =~ s/f2=//;    
    $f2 = $argv;
    }  
  }

$filein =~ tr/A-Z/a-z/;
$fileout = $filein;
$fasta = $filein;
$seg = $filein;
$gene = $filein;
$gene =~ s/_.*//g;

print "\@ARGV: \$filein = '$filein', \$remove_mono = '$remove_mono', \$seq_children = '$seq_children', \$omit = '$omit' \$f2=$f2\n";

&get_starting_seqs;

# Squishing - Iterative process to remove monomorphic positions, identify SEG groups and remove child-alleles

$no_alleles[0] = @allele_index;

# Build List of Alleles to include in iteration
  
@allele_index = ();  
for $allele ( keys %alleles )
  {
  if ( $alleles{$allele}{'analyse'} eq 1 )
    {
    push(@allele_index, $allele);    
    }  
  }
$no_alleles = @allele_index;  
$number_alleles = @allele_index;  
print "AI: $number_alleles (pre-SNP filter)\n";
 
# Identify alleles differing by a single SNP 
 
if ( $remove_snps == 0 )
  {
  #@allele_index = sort_new_nomenc::sorting_nomenclature(@allele_index);
  for ( $a=0; $a<@allele_index; ++$a )
    {
    $allele_key_a = $allele_index[$a];
    for ( $b=$a+1; $b<@allele_index; ++$b )
      {
      $allele_key_b = $allele_index[$b];
      
      $num_mismatch = hd($alleles{$allele_key_a}{'sequence'}, $alleles{$allele_key_b}{'sequence'});
      if ( $num_mismatch < 2 )
        {
        $alleles{$allele_key_b}{'analyse'} = 0;  
        }
      }    
    }
  }
  
@allele_index = ();  
for $allele ( keys %alleles )
  {
  if ( $alleles{$allele}{'analyse'} eq 1 )
    {
    push(@allele_index, $allele);    
    }  
  }

$number_alleles = @allele_index;  
print "AI: $number_alleles (post snp filter)\n";

print "$allele_index[0] $alleles{$allele_index[0]}{'sequence'}\n";
@base_keys = ('A','C','G','T','D','I');
$no_alleles = @allele_index;
$bases_per_position = ();

open(RDAT, ">r_in_ex2.dat");
print RDAT "Position\tBases\tF2\n";
for ( $bp=0; $bp<270; ++$bp )
  {
  %base_count = ();
  $position = $bp+$pos_offset;
  
  foreach $allele (@allele_index)
    {
    $base = substr($alleles{$allele}{'sequence'}, $bp, 1);    
    ++$base_count{$base};
    }
  
  @bases = keys(%base_count);
  $base_count = @bases;
  
  $bases_per_position{$base_count}++;
  
  @f2_values = ();
  foreach $base_key (@base_keys)
    {
    $f2_value = $base_count{$base_key}/$no_alleles*100;
    push(@f2_values, $f2_value);
    }
  
  @f2_values = sort { $b <=> $a } @f2_values;
  
  print "$position\t$f2_values[1]\n";
  
  --$base_count;
  $base_count =~ s/0/-1/;
  
  if ( $position%30 eq 0)
    {
    print RDAT "$position\t$base_count\t$f2_values[1]\n";
    }
  else
    {
    print RDAT "\t$base_count\t$f2_values[1]\n";
    }
  }

close(RDAT);

print "Number of Nucleotide Bases seen per position E2\n";

foreach $count ( sort ( keys (%bases_per_position) ) )
  {
  print "$count\t $bases_per_position{$count}\n";     
  }


#exit;

open(RDAT, ">r_in_ex3.dat");
print RDAT "Position\tBases\tF2\n";
for ( $bp=270; $bp<546; ++$bp )
  {
  %base_count = ();
  $position = $bp+$pos_offset;
  
  foreach $allele (@allele_index)
    {
    $base = substr($alleles{$allele}{'sequence'}, $bp, 1);    
    ++$base_count{$base};
    }
  
  @bases = keys(%base_count);
  $base_count = @bases;
  
  $bases_per_position{$base_count}++;
  
  @f2_values = ();
  foreach $base_key (@base_keys)
    {
    $f2_value = $base_count{$base_key}/$no_alleles*100;
    push(@f2_values, $f2_value);
    }
  
  @f2_values = sort { $b <=> $a } @f2_values;
  
  print "$position\t$f2_values[1]\n";
  
  --$base_count;
  $base_count =~ s/0/-1/;
  
  if ( $position%30 eq 0)
    {
    print RDAT "$position\t$base_count\t$f2_values[1]\n";
    }
  else
    {
    print RDAT "\t$base_count\t$f2_values[1]\n";
    }
  }
close(RDAT);

$title_gene = $gene;
$title_gene =~ tr/a-z/A-Z/;
if ( $remove_snps == 0 )
  {
  $title_desc = 'Motifs';
  $file_label = '_motifs';
  }
else
  {
  $title_desc = 'SNPs';  
  $file_label = '_snps';
  }
 
print "USE FILE LABEL: '$file_label'\n"; 
  
use Statistics::R;

print "Building $gene\_ex2$file_label\.png\n";

$R = Statistics::R->new( shared => 1);
$R->run(qq`snp_data <- read.table("r_in_ex2.dat", header=T, sep="\t")`);
$R->run(qq`png("$gene\_ex2$file_label\.png", width=1600)`);
$R->run(qq`barplot(snp_data\$Bases, names.arg = snp_data\$Position, lwd=2, ylim=c(-1,4), axes=FALSE, beside=TRUE, xpd=FALSE, xaxs="i", las=1, cex.names=2, cex.label=2, cey.names=2, cey.label=2, col.axis="black", col="black")`);
$R->run(qq`axis(side=2,at=c(-1,0,1,2,3,4), labels=c('1','','2','3','4','5'), lwd=2, cex.axis=2, las=1)`);
$R->run(qq`dev.off()`);
$R->stop();

print "Building $gene\_ex2\_f2_freq\.png\n";

$R=Statistics::R->new( shared => 1);
$R->run(qq`snp_data <- read.table("r_in_ex2.dat", header=T, sep="\t")`);
$R->run(qq`png("$gene\_ex2\_f2_freq\.png", width=1600)`);
$R->run(qq`barplot(snp_data\$F2, names.arg=snp_data\$Position, ylim=c(0,50), axes=FALSE, beside=TRUE, xpd=FALSE, xaxs="i", las=1, lwd=2, cex.names=2, cex.label=2, cex.axis=2, col.axis="black", col="black")`);
$R->run(qq`axis(side=2, lwd=2, cex.axis=2, las=1)`);
$R->run(qq`dev.off()`);
$R->stop();

print "Building $gene\_ex3$file_label\.png\n";

$R=Statistics::R->new( shared => 1);
$R->run(qq`snp_data <- read.table("r_in_ex3.dat", header=T, sep="\t")`);
$R->run(qq`png("$gene\_ex3$file_label\.png", width=1600)`);
$R->run(qq`barplot(snp_data\$Bases, names.arg=snp_data\$Position, lwd=2, ylim=c(-1,4), axes=FALSE, beside=TRUE, xpd=FALSE, xaxs="i", las=1, cex.names=2, cex.axis=2, cey.names=2, cey.axis=2, col.axis="black", col="black")`);
$R->run(qq`axis(side=2,at=c(-1,0,1,2,3,4),labels=c('1','','2','3','4','5'), lwd=2, cex.axis=2, las=1)`);
$R->run(qq`dev.off()`);
$R->stop();

print "Building $gene\_ex3\_f2_freq\.png\n";

$R=Statistics::R->new( shared => 1);
$R->run(qq`snp_data <- read.table("r_in_ex3.dat", header=T, sep="\t")`);
$R->run(qq`png("$gene\_ex3\_f2_freq\.png", width=1600)`);
$R->run(qq`barplot(snp_data\$F2, names.arg=snp_data\$Position, ylim=c(0,50), axes=FALSE, beside=TRUE, xpd=FALSE, xaxs="i", las=1, lwd=2, cex.names=2, cex.label=2, cey.names=2, cey.label=2, col.axis="black", col="black")`);
$R->run(qq`axis(side=2, lwd=2, cex.axis=2, las=1)`);
$R->run(qq`dev.off()`);
$R->stop();

print "Number of Nucleotide Bases seen per position\n";

foreach $count ( sort ( keys (%bases_per_position) ) )
  {
  print "$count\t $bases_per_position{$count}\n";     
  }


exit;

sub hd
  {
  return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
  }

sub define_substitutions  
  {
  $allele{'HLA-A_31:57'}{'substitute'} = 'HLA-A_31:01:02';
  
  $allele{'HLA-B_07:08'}{'substitute'}  = 'HLA-B_07:02:01';
  $allele{'HLA-B_08:03'}{'substitute'}  = 'HLA-B_08:01:01';
  $allele{'HLA-B_45:05'}{'substitute'}  = 'HLA-B_45:01';
  $allele{'HLA-B_51:10'}{'substitute'}  = 'HLA-B_51:01:01';
  $allele{'HLA-B_35:167'}{'substitute'} = 'HLA-B_35:168';
  
  $allele{'HLA-C_07:65'}{'substitute'}  = 'HLA-C_07:60';
  $allele{'HLA-C_12:60'}{'substitute'}  = 'HLA-C_12:02:01';
  }
  
sub get_starting_seqs
  {
  print "Using sequences from $filein\n";
    
  open(LIB, $filein);
  while(<LIB>)
    {  
    chomp();
    if ( $_ =~ /,/ )
      {
      next;    
      }
    $_ =~ s/[\n\r]//;  
    $_ =~ s/\t[NY]//;  
    $_ =~ s/\|//;  
    ($allele,$seq) = split(/\t/, $_);
  
    $refseq = $seq;
    last;
    }
  close(LIB);
    
  open(LIB, "$filein");
  while(<LIB>)
    {
    chomp();
    if ( $_ =~ /,/ )
      {
      next;    
      }
    $_ =~ s/[\n\r]//;  
    $_ =~ s/\t[NY]//;  
    $_ =~ s/\|//;  

    ($allele,$seq) = split(/\t/, $_);
    $allele =~ s/\*/\_/;
    $allele =~ s/>//;
  
    if ( ( $omit =~ /nhp/ || $omit =~ /NHP/ ) && $allele !~ /^HLA/ )
      {
      next;  
      }
    #if ( ( $omit =~ /nhp/ || $omit =~ /NHP/ ) && $allele =~ /^HLA/ )
    #  {
    #  $test340 = `grep ^$allele\$ hla_nom.txt`;
    #  chomp($test340);
    #  if ( $test340 ne $allele )
    #    {
    #    next;
    #    }
    #  }  
    if ( ( $omit =~ /hla/ || $omit =~ /HLA/ ) && $allele =~ /^HLA/ )
      {
      next;  
      }
    
    $alleles{$allele}{'sequence'} = $seq;
    $alleles{$allele}{'analyse'}  = 1;    
    }
  close(LIB);
    
  @allele_index = keys(%alleles); 
  $allele_index = \@allele_index;
  
  # Recode InDels to either I or D. I = indel plus following base all coded as a single I, 
  
  for ( $bp=0; $bp<length($refseq); ++$bp )
    {
    $refbase = substr($refseq, $bp, 1);  
    if ( $refbase eq '.' )
      {
      substr($refseq, $bp, 1, '!');
      for ( $a=0; $a<@allele_index; ++$a )
        {
        $allele_key = $allele_index[$a];
        $seqbase =  substr($alleles{$allele_key}{'sequence'}, $bp, 1); 
        if ( $seqbase eq '.' )  
          {
          substr($alleles{$allele_key}{'sequence'}, $bp, 1, '!');
          }
        elsif( $seqbase =~ /[ACTGX]/ )  
          {
          substr($alleles{$allele_key}{'sequence'}, $bp, 1, 'I');
          }
        }
      }
    }
  
  for $allele ( keys %alleles )
    {
    $alleles{$allele}{'sequence'} =~ s/I{1,}[ACGTX]/I/g;  
    $alleles{$allele}{'sequence'} =~ s/\./D/g;  
    $alleles{$allele}{'sequence'} =~ s/!//g;
    
    if ( $alleles{$allele}{'sequence'} =~ /[ID]/ )
      {
      #print "Omitting $allele\n";
      }
    }
  #exit;
  $number_alleles = @allele_index;  
  print "AI: $number_alleles (initial load)\n";
  
  # Substitute unsequenced regions
  
  &define_substitutions;
  
  for $allele ( keys %alleles )
    {
    if ( $alleles{$allele}{'sequence'} =~ /\*/ )
      {
      print "WARNING! Unsequenced region remains in $allele[$a]\n";     
      }
    }
  }

sub hd
  {
  return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
  }
 


