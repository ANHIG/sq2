#!/usr/bin/perl
#
# A perl script to generate variation dotplots as described in;
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

use sort_new_nomenc;
use remove_monomorphic;

$remove_mono  = 0;
$omit         = 0;
$seg_children = 1;
$pos_offset   = 74;

foreach $argv (@ARGV) 
  {
  if ( $argv =~ /^i/ )
    {
    $argv =~ s/i=//;    
    $filein = $argv;
    }
  elsif ( $argv =~ /^r/ )
    {
    $argv =~ s/r=//;    
    $remove_mono = $argv;
    }
  elsif ( $argv =~ /^c/ )
    {
    $argv =~ s/c=//;    
    $children = $argv;
    }  
  elsif ( $argv =~ /^h/ )
    {
    $argv =~ s/h=//;    
    $histo_group = $argv;
    }  
  elsif ( $argv =~ /^o/ )
    {
    $argv =~ s/o=//;    
    $omit = $argv;
    }  
  }

$filein =~ tr/A-Z/a-z/;
$fileout = $filein;
$fasta = $filein;
$seg = $filein;

print "\@ARGV: \$filein = '$filein', \$remove_mono = '$remove_mono', \$seq_children = '$seq_children', \$omit = '$omit'\n";

&get_starting_seqs;

# Squishing - Iterative process to remove monomorphic positions, identify SEG groups and remove child-alleles

$no_alleles[0] = @allele_index;

print "Iteration 0: $no_alleles[0] (Alleles loaded)\n";

for ( $i=1; $i<99; ++$i )
  {
  if ( $no_alleles[$i-1] eq $no_alleles[$i] )
    {
    return;    
    }
    
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
  
  # Remove Identical Sequences from current iteration
  
  @allele_index = sort_new_nomenc::sorting_nomenclature(@allele_index);
  for ( $a=0; $a<@allele_index; ++$a )
    {
    $allele_key_a = $allele_index[$a];
    if ( $alleles{$allele_key_a}{'analyse'} eq 0 )
      {
      next;
      }
    for ( $b=$a+1; $b<@allele_index; ++$b )
      {
      $allele_key_b = $allele_index[$b];
      if ( $alleles{$allele_key_a}{'sequence'} eq $alleles{$allele_key_b}{'sequence'} )
        {
        $alleles{$allele_key_b}{'analyse'} = 0;
        $iteration{$i}{$allele_key_b} = 'Removed, identical to ' . $allele_key_a;
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
  $no_alleles = @allele_index;
  
  print "\tDistinct Alleles: $no_alleles\n";
 
  # Identify alleles differing by a single SNP 
 
  @allele_index = sort_new_nomenc::sorting_nomenclature(@allele_index);
  for ( $a=0; $a<@allele_index; ++$a )
    {
    $allele_key_a = $allele_index[$a];
    for ( $b=0; $b<@allele_index; ++$b )
      {
      $allele_key_b = $allele_index[$b];
      
      $num_mismatch = hd($alleles{$allele_key_a}{'sequence'}, $alleles{$allele_key_b}{'sequence'});

      if ( $num_mismatch eq 1 && $seg_groups{$allele_key_a}{'children'}[0] !~ /[0-9]/ )
        {
        $seg_groups{$allele_key_a}{'children'}[0] = $allele_key_b;
        $seg_groups{$allele_key_a}{'size'} = 1;
        }    
      elsif ( $num_mismatch eq 1 )
        {
        push @{$seg_groups{$allele_key_a}{'children'}}, $allele_key_b;
        ++$seg_groups{$allele_key_a}{'size'};
        }
      }  
    }

  # Analyse SEG groups, stage 1 sort alleles within the SEG, then identify SEGs of identical composition and remove duplicates
  print "Sorting SEG groups members\n";

  @allele_index = ();  
  
  for $allele (sort(keys %alleles))
    {
    @f1child = @{$seg_groups{$allele}{'children'}};
    $no_f1 = @f1child;
    foreach $f1child ( @f1child )
      {
      #print "$allele,$f1child,$no_f1\n";
      }
    }
    
  print "Point Mutants include:\n";
    
  for $allele (sort(keys %alleles))
    {
    @f1child = @{$seg_groups{$allele}{'children'}};
    $no_f1 = @f1child;
    
    print "@f1child\n";
    
    if ( $no_f1 < 1 )
      {
      $alleles{$allele}{'analyse'} = 0;
      }
    
    foreach $f1child ( @f1child )
      {
      if ( !defined($seg_groups{$f1child}) )
        {
        next;
        }
      @f2child = @{$seg_groups{$f1child}{'children'}};
      $no_f2 = @f2child;
    
      if ( $no_f2 < $no_f1 )
        {
        #print "\t$f1child is child allele\n";
        $alleles{$f1child}{'analyse'} = 0;      
        }
      if ( ( $no_f1 == 1 && $no_f2 == 1 && $alleles{$allele}{'analyse'} ne 0 ) )
        {
        #print "Recipricol pairing $allele ( $f1child ($alleles{$f1child}{'analyse'}) <> $f2child[0] ($alleles{$f2child[0]}{'analyse'}) )\n";
        $alleles{$f1child}{'analyse'} = 0;
        #$alleles{$f2child[0]}{'analyse'} = 0;
        }
      }
    
    }
  exit;
   
  #%seq = ();
  # 
  
  $point_mutant_total = 0;
  
  @allele_index = ();  
  for $allele (sort ( keys %alleles ) )
    {
    if ( $alleles{$allele}{'analyse'} eq 1 )
      {
      ++$count;
    
      @f1child = @{$seg_groups{$allele}{'children'}};
      $no_f1 = @f1child;
      $point_mutant_total = $point_mutant_total+$no_f1;
      
      print "$count $allele ($no_f1) [@f1child]\n"; 
      push(@allele_index, $allele);
  #    $seq{$allele} = $alleles{$allele}{'sequence'};
  
      foreach $f1child (@f1child)
        {
        #print "$allele,$f1child\n";
        }
       
      }
    }
  $no_alleles = @allele_index;    

print "@allele_index\n";
print "Point Mutant Total = $point_mutant_total\n";
  
  print "Generated Lineages for SEGs, $no_alleles alleles left for analysis\n";
 
  # SEG CHECK
  
  for ( $a=0; $a<@allele_index; ++$a )
    {
    $allele_key_a = $allele_index[$a];
    for ( $b=0; $b<@allele_index; ++$b )
      {
      $allele_key_b = $allele_index[$b];
      
      $num_mismatch = hd($alleles{$allele_key_a}{'sequence'}, $alleles{$allele_key_b}{'sequence'});

      if ( $num_mismatch eq 1 && $alleles{$allele_key_a}{'analyse'} ne 0 )
        {
        print "WARNING SEG STILL FOUND FOR $allele_key_a vs $allele_key_b removing $allele_key_b SEG\n"; 
        $alleles{$allele_key_b}{'analyse'} = 0;  
        #exit;
        }
      }  
    }
    
  @allele_index = ();  
  for $allele (sort ( keys %alleles ) )
    {
    if ( $alleles{$allele}{'analyse'} eq 1 )
      {
      ++$count;
      push(@allele_index, $allele);
      }
    }
  $no_alleles = @allele_index;    

  print "@allele_index\n";  
    
  print "\nMoving to motif analysis\n\n";
  
  # Remove monomorphic positions
  
  for ( $a=0; $a<@allele_index; ++$a )
    {
    $motif_seq{$allele_index[$a]} = $alleles{$allele_index[$a]}{'sequence'};
    }
    
  for ( $bp=0; $bp<length($motif_seq{$allele_index[0]}); ++$bp )
    {
    $refbase = substr($motif_seq{$allele_index[0]}, $bp, 1);  
    $mono = 1;
    for ( $t=1; $t<@allele_index; ++$t )
      {
      $seqbase = substr($motif_seq{$allele_index[$t]}, $bp, 1);  
      if ( $seqbase ne $refbase )
        {
        $mono = 0;
        push(@varpos, $bp+74);
        last;
        }
      }
    if ( $mono eq 1 )
      {
      for ( $t=0; $t<@allele_index; ++$t )
        {
        substr($motif_seq{$allele_index[$t]}, $bp, 1, '!');
        }
      }
    }
  
  for ( $t=0; $t<@allele_index; ++$t )
    {
    $motif_seq{$allele_index[$t]} =~ s/!//g;
    }


  print "Variable positions at @varpos\n";

  use Spreadsheet::WriteExcel;
  my $workbook = Spreadsheet::WriteExcel->new('core_sq2_v3.xls');
  $worksheet1 = $workbook->add_worksheet('Potential Motifs');
  
  $seqlength = length($motif_seq{$allele_index[0]});
  
  $row = 0;
  foreach $allele (@allele_index)
      {
      ++$row;
      $worksheet1->write($row, 1, $allele);
      }
  
  $col=1;   
  for ( $bp=0; $bp<$seqlength; ++$bp )
    {
    %base_count = ();
    foreach $allele (@allele_index)
      {
      $base = substr($motif_seq{$allele}, $bp, 1);
      $base_count{$base}++;
      } 
    @base_count = keys(%base_count);

    if ( @base_count > 1 )
      {
      #print "$bp @base_count\n";
      ++$col;
      $worksheet1->write(0, $col, $varpos[$bp]);

      $row = 0;
      foreach $allele (@allele_index)
      {
      ++$row;
      
      $base = substr($motif_seq{$allele}, $bp, 1);
      if ( $base eq 'A' )
          {
          $format = $workbook->add_format();
          $format->set_format_properties(bg_color => 'green');
          $worksheet1->write($row, $col, $base, $format);
          }
        elsif ( $base eq 'C' )
          {
          $format = $workbook->add_format();
          $format->set_format_properties(bg_color => 'blue');
          $worksheet1->write($row, $col, $base, $format);
          }
        elsif ( $base eq 'G' )
          {
          $format = $workbook->add_format();
          $format->set_format_properties(bg_color => 'black', color => 'white');
          $worksheet1->write($row, $col, $base, $format);
          }
        elsif ( $base eq 'T' )
          {
          $format = $workbook->add_format();
          $format->set_format_properties(bg_color => 'red');
          $worksheet1->write($row, $col, $base, $format);
          }
       }
     
    }
    
  }
    
  # Look at shared motifs, using a rolling window of 10bps every 5 bps. Score each block and group identical units. Then review overlap.
  
 tie @shared_motifs, 'Array::Unique';
 tie @meg_motifs, 'Array::Unique';
 @shared_motifs = ();
 @shared_motifs = sort( @shared_motifs );

 foreach $allele_a ( @allele_index )
    {
    ++$a;
    $min_diffs = 99;  
    for ( $b=0; $b<@allele_index; ++$b )
      {
      $allele_b = $allele_index[$b];  
      if ( substr($allele_a, 0, 8) eq substr($allele_b, 0, 8) )
        {
        next;
        }
        
      $diff = hd( $motif_seq{$allele_a}, $motif_seq{$allele_b} );
      
      #print "$allele_a vs $allele_b = $diff\n";
      if ( $diff < $min_diffs )
        {
        $min_diffs = $diff;
        $min_allele = $allele_b;
        }
      }
    $shared_motif = '';
    $bp=0;
    
    #print "\tidentify diff of $min_difffs between $allele_a and $min_allele\n";
    
    @motpos = ();
    for ( $bp=0; $bp<length($motif_seq{$allele_a}); ++$bp )
      {
      if ( substr($motif_seq{$allele_a}, $bp, 1) ne substr($motif_seq{$min_allele}, $bp, 1) )
        {
        push(@motpos, $bp);
        }  
      }
      
    $vpst = $motpos[0];
    $vpe  = $motpos[-1];
    $vlen = length($motif_seq{$allele_a});
    print "DEBUG $allele_a vs $min_allele (Based on $motpos[0] to $motpos[-1] of $vlen)\n";
    
    $conserved_sequence1 = $motif_seq{$allele_a};
    $blanklength = $vlen - $vpst;
    $blank = ' ' x $blanklength;
    $conserved_sequence1 =~ s/(^[ACTG]{$vpst})([ACTG].*)/$1$blank/;
    
    $conserved_sequence2 = $motif_seq{$allele_a};
    $blanklength = $vpe;
    $blank = ' ' x $blanklength;
    $conserved_sequence2 =~ s/(^[ACTG]{$blanklength})([ACTG].*$)/$blank$2/;
    
    $recombinant_sequence1 = $motif_seq{$allele_a};
    $blanklength1 = $vpst;
    $blanklength2 = $vlen-$vpe;
    $blank1 = ' ' x $blanklength1;
    $blank2 = ' ' x $blanklength2;
    $recombinant_sequence1 =~ s/(^[ACTG]{$blanklength1})([ACTG].*)([ACTG]{$blanklength2}$)/$blank1$2$blank2/;

    $recombinant_sequence2 = $motif_seq{$min_allele};
    $blanklength1 = $vpst;
    $blanklength2 = $vlen-$vpe;
    $blank1 = ' ' x $blanklength1;
    $blank2 = ' ' x $blanklength2;
    $recombinant_sequence2 =~ s/(^[ACTG]{$blanklength1})([ACTG].*)([ACTG]{$blanklength2}$)/$blank1$2$blank2/;   
        
    push(@meg_motifs, $conserved_sequence1);
    push(@meg_motifs, $recombinant_sequence1);
    push(@meg_motifs, $recombinant_sequence2);
    push(@meg_motifs, $conserved_sequence2);
   }
   
   
 print "Shared Motif Index\n";
 $view_varpos = join("\t", @varpos);
 print "Variable positions\n";
 print "\n";
 foreach $varpos ( @varpos )
   {
   $varposf = sprintf("%03d", $varpos);
   $varposf =~ s/^0/ /;
   $varposp = substr($varposf, 0, 1);
   print "$varposp";
   }
 print "\n";
 foreach $varpos ( @varpos )
   {
   $varposf = sprintf("%03d", $varpos);
   $varposf =~ s/^0/ /;
   $varposp = substr($varposf, 1, 1);
   print "$varposp";
   }
 print "\n";
 foreach $varpos ( @varpos )
   {
   $varposf = sprintf("%03d", $varpos);
   $varposf =~ s/^0/ /;
   $varposp = substr($varposf, 2, 1);
   print "$varposp";
   }
 print "\n";

   @meg_motifs = sort ( @meg_motifs);
   foreach $meg ( @meg_motifs )
     {
     if ( $count = $meg =~ tr/[ACGT]// > 1 )
       {
       print "MEG '$meg'\n";
       @msa_lookup = $meg =~ /(^ {0,})([ACTG]*)( {0,}$)/;
       $msa_lookup_start = length($msa_lookup[0]);
       $msa_lookup_end   = length($msa_lookup[2]);
   	   
       %meg_hits = ();
       for $allele_c (@allele_index )
         {
         if ( substr($motif_seq{$allele_c}, $msa_lookup_start, ) =~ /$msa_lookup[1]/ )
           {
           $allele_c =~ s/:.*//;
           ++$meg_hits{$allele_c};
           }
         #exit;
         }
       foreach $key ( keys %meg_hits )
         {
         print "$key ($meg_hits{$key}) ";
         }
       print "\n";
       }
     }

 @meg_motifs = sort(@meg_motifs);   
 
 print "Shared Motif Index\n";
 $view_varpos = join("\t", @varpos);
 print "Variable positions\n";
 print "\n\t";
 foreach $varpos ( @varpos )
   {
   $varposf = sprintf("%03d", $varpos);
   $varposf =~ s/^0/ /;
   $varposp = substr($varposf, 0, 1);
   print "$varposp";
   }
 print "\n\t";
 foreach $varpos ( @varpos )
   {
   $varposf = sprintf("%03d", $varpos);
   $varposf =~ s/^0/ /;
   $varposp = substr($varposf, 1, 1);
   print "$varposp";
   }
 print "\n\t";
 foreach $varpos ( @varpos )
   {
   $varposf = sprintf("%03d", $varpos);
   $varposf =~ s/^0/ /;
   $varposp = substr($varposf, 2, 1);
   print "$varposp";
   }
 print "\n";


 foreach $meg ( @meg_motifs )
   {
   print "$meg\n";
   }
  
  $a = 0;
  foreach $shared_motif ( @shared_motifs )
    {
    ++$a;  
    for ( $b=$a+1; $b<@shared_motifs; ++$b )
      {
      $shared_motif_b = $shared_motifs[$b];
      
      if (  ( $shared_motif =~ $shared_motif_b || $shared_motif_b =~ $shared_motif ) &&
          join("", @{$motif_matches{$shared_motif}{'alleles'}}) eq join("",@{$motif_matches{$shared_motif_b}{'alleles'}}) )
        {
        #print "Comparing $shared_motif ($motif_matches{$shared_motif}{'size'}) vs $shared_motif_b ($motif_matches{$shared_motif_b}{'size'})\n";
        if ( $motif_matches{$shared_motif}{'size'} < $motif_matches{$shared_motif_b}{'size'} )
          {
          delete($motif_matches{$shared_motif});  
          print "\tDeleting $shared_motif\n";
          }
        elsif ( $motif_matches{$shared_motif}{'size'} > $motif_matches{$shared_motif_b}{'size'} )
          {
          delete($motif_matches{$shared_motif_b});  
          print "\tDeleting $shared_motif_b\n";
          }
        }
      }
    }
  
  @shared_motifs = sort ( keys %motif_matches );
  $no_shared_motifs = @shared_motifs;
  
  # Looking for core postions within Motifs, some based are included in multiple motifs. What are these core motifs and can we use these to rank the motifs.
  
  foreach $shared_motif ( sort @shared_motifs )
    {
    #print "$shared_motif\t$motif_matches{$shared_motif}{'Child'}\t$motif_matches{$shared_motif}{'Parent'}\t$motif_matches{$shared_motif}{'size'}\t$motif_matches{$shared_motif}{'length'}\t@{$motif_matches{$shared_motif}{'alleles'}}\n";  
    ++$parental{$motif_matches{$shared_motif}{'Parent'}};
    }  
  
  foreach $key ( keys %parental )
    {
    print "KEY: '$key'\tVALUE: '$parental{$key}'\n";  
    }
  
  
  $rank = 0;
  foreach $parent ( sort { $parental{$b} <=> $parental{$a} } ( keys %parental ) )
    {
    ++$rank;
    $rank{$parent} = $rank;
    print "$parent Rank: $rank{$parent}\n";
    }
    
  print "\n\n";  
    
  foreach $parent ( sort { $parental{$b} <=> $parental{$a} } ( keys %parental ) )
    {
    print "$parent Stored Rank: $rank{$parent}\n";
    
    foreach $shared_motif ( @shared_motifs )
      {
      if ( $motif_matches{$shared_motif}{'Child'} eq $parent )
        {
        $p = $motif_matches{$shared_motif}{'Parent'};
        $c = $motif_matches{$shared_motif}{'Child'};
        #print "\tGRANDCHILD?: $shared_motif\tP: $p\tC: $c\n";
        #print "\t$p $rank{$p} > $c $rank{$c}\n";
        if ( $rank{$p} > $rank{$c} )
          {
          $omit_core{$p} = 'Y';  # Problems with the parental hash??????
          }
        print "\tOMIT: $p $omit_core{$p}\n";
        }
      }
    }
  
  print "\n\n";
  
  exit;
  
  foreach $parent  ( keys %parental ) 
    {
    if ( $omit_core{$parent} !~ 'Y' )
      {
      print "CORE ALLELE: '$parent'\n"; 
      }
    }
    
   
    
    
    
  foreach $core_position ( sort { $core_positions{$b} <=> $core_positions{$a} } ( keys %core_positions ) )
    {
    print "Looking for core positions in the motif: $core_position ($core_positions{$core_position})\n";
    }
    
  
    
  foreach $shared_motif ( sort { $motif_matches{$b}{'AVG Score'} <=> $motif_matches{$a}{'AVG Score'}} @shared_motifs )
    {
    print "$shared_motif\t$motif_matches{$shared_motif}{'size'}\t$motif_matches{$shared_motif}{'Motif Score'}\t$motif_matches{$shared_motif}{'length'}\t$motif_matches{$shared_motif}{'AVG Score'}\t@{$motif_matches{$shared_motif}{'alleles'}}\n";  
    }  
  
  print "No. of Shared Motifs: $no_shared_motifs\n";
  
  # Analyse SEGs and remaining alleles for SEGs based on motifs as single SNPs
  # Foreach motif, append sequence with an M if motif present, and N if not. Remove all bases of motif from the sequence before calculating HD.
  
  @allele_index = sort_new_nomenc::sorting_nomenclature(@allele_index);
  
  foreach $shared_motif ( @shared_motifs )
    {
    #print "Looking for SEGs with Motifs for $shared_motif (@{$motif_matches{$shared_motif}{'alleles'}})\n";
    @search = split(/:/, $shared_motif);
      
    $a=0;
    foreach $allele_key_a ( @allele_index )
      {
      $seq_a = $alleles{$allele_key_a}{'sequence'};
      foreach $snp ( @search )
        {
        $position = $snp;
        substr($seq_a, $position-$pos_offset, 1, " ");
        }

      if ( grep ( /^$allele_key_a$/, @{$motif_matches{$shared_motif}{'alleles'}} ) )
        {
        $seq_a =~ s/ /M/;
        }
      else
        {
        $seq_a =~ s/ /N/;    
        }
      $seq_a =~ s/ //g;  

      for ( $b=$a+1; $b<@allele_index; ++$b )
        {
        $allele_key_b = $allele_index[$b];
        $seq_b = $alleles{$allele_key_b}{'sequence'};
        foreach $snp ( @search )
          {
          $position = $snp;
          substr($seq_b, $position-$pos_offset, 1, " ");
          }

        if ( grep ( /^$allele_key_b$/, @{$motif_matches{$shared_motif}{'alleles'}} ) )
          {
          $seq_b =~ s/ /M/;
          }
        else
          {
          $seq_b =~ s/ /N/;    
          }
        $seq_b =~ s/ //g;  
        $motif_hd =  hd($seq_a, $seq_b); 
         
        if ( $motif_hd eq 1 && $seq_a =~ /M/ && $seq_b =~ /M/ )
          {
          print "$shared_motif (@{$motif_matches{$shared_motif}{'alleles'}})\nA: $allele_key_a:\t$seq_a\nB: $allele_key_b:\t$seq_b\n\n";  
          #exit;
          }     
        }
      }  
    }
  
  
  #exit;
  
  # Sorted by size
  
  for $shared_motif ( sort { $motif_matches{$b}{'size'} <=> $motif_matches{$a}{'size'} } keys %motif_matches  )
    {
    print "Alleles for $shared_motif: $motif_matches{$shared_motif}{'length'} bps: $motif_matches{$shared_motif}{'size'}/$no_alleles[$i] :\n@{$motif_matches{$shared_motif}{'alleles'}}\n";
    }  
  
  # Sorted by length
    
  for $shared_motif ( sort { $motif_matches{$b}{'length'} <=> $motif_matches{$a}{'length'} } keys %motif_matches  )
    {
    #print "Alleles for $shared_motif: $motif_matches{$shared_motif}{'length'} bps: $motif_matches{$shared_motif}{'size'}/$no_alleles[$i] :\n@{$motif_matches{$shared_motif}{'alleles'}}\n";
    }  
    
  # Unsorted
    
  for $shared_motif ( keys %motif_matches  )
    {
    #print "$shared_motif\t$motif_matches{$shared_motif}{'length'}\t$motif_matches{$shared_motif}{'size'}\t@{$motif_matches{$shared_motif}{'alleles'}}\n";
    }  
  
  exit;
  }

exit;

sub hd
  {
  return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
  }

sub define_substitutions  
  {
  $alleles{'HLA-A_31:57'}{'substitute'} = 'HLA-A_31:01:02';
  
  $alleles{'HLA-B_07:08'}{'substitute'}  = 'HLA-B_07:02:01';
  $alleles{'HLA-B_08:03'}{'substitute'}  = 'HLA-B_08:01:01';
  $alleles{'HLA-B_45:05'}{'substitute'}  = 'HLA-B_45:01';
  $alleles{'HLA-B_51:10'}{'substitute'}  = 'HLA-B_51:01:01';
  $alleles{'HLA-B_35:167'}{'substitute'} = 'HLA-B_35:168';
  
  $alleles{'HLA-C_07:65'}{'substitute'}  = 'HLA-C_07:60';
  $alleles{'HLA-C_12:60'}{'substitute'}  = 'HLA-C_12:02:01';
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
      #print "SKIPPING $allele - NHP in HLA DATASET\n";  
      next;  
      }
    if ( ( $omit =~ /hla/ || $omit =~ /HLA/ ) && $allele =~ /^HLA/ )
      {
      #print "SKIPPING $allele - HLA in NHP DATASET\n";  
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
    
    #print "$allele ($alleles{$allele}{'analyse'})-> $alleles{$allele}{'sequence'}\n";
    }
    
  # Substitute unsequenced regions
  
  &define_substitutions;
  
  for $allele ( keys %alleles )
    {
    if ( $alleles{$allele}{'sequence'} =~ /\*/ )
      {
      print "WARNING! Unsequenced region remains in $allele ...";
      $sub = $alleles{$allele}{'substitute'};

      for ( $bp=0; $bp<length($alleles{$allele}{'sequence'}); ++$bp )
        {
        $par = substr($alleles{$allele}{'sequence'}, $bp, 1);
        $ext = substr($alleles{$sub}{'sequence'}, $bp, 1);
        
        if ( $par =~ /\*/ )
          {
          substr($alleles{$allele}{'sequence'}, $bp, 1, $ext); 
          }   
        }
      print " padded using $alleles{$allele}{'substitute'}\n";
      }
    }
  }

sub con_base
  {
  @counts = ();  
  @counts = reverse(sort { $a <=> $b } (values(%count)));
  
  #print " -> @counts\n";
  foreach $key (keys %count)
    {
    if ( $count{$key} eq $counts[0] )
      {
      $conbase = $key;  
      }
    } 
  $consensus = join(':', $consensus, $conbase);
  %count = ();  
  $count{'A'}=$count{'C'}=$count{'G'}=$count{'T'} = 0;
  }

sub hd
  {
  return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
  }
 


