#!/usr/bin/perl
#
# Module to remove monomorphic positions from aligned sequences
#-------------------------------------------------------------------------------

package remove_monomorphic;

sub processing ( % )
  {
  use sort_new_nomenc;

  $ref = shift;
  
  @to_process = keys(%{$ref});
  @to_process = sort_new_nomenc::sorting_nomenclature(@to_process);
  pop(@to_process);
  
  %vals = %{$ref};  
  $positions = delete($vals{'POSITIONS'});
  @positions = split(/,/, $positions);
  
  for ( $bp=0; $bp<length($vals{$to_process[0]}); ++$bp )
    {
    $refbase = substr($vals{$to_process[0]}, $bp, 1);  
    $mono = 1;
    for ( $t=1; $t<@to_process; ++$t )
      {
      $seqbase = substr($vals{$to_process[$t]}, $bp, 1);  
      if ( $seqbase ne $refbase )
        {
        $mono = 0;
        last;
        }
      }
    if ( $mono eq 1 )
      {
      $positions[$bp] = '!';  
      for ( $t=0; $t<@to_process; ++$t )
        {
        substr($vals{$to_process[$t]}, $bp, 1, '!');
        }
      }
    }
  
  for ( $t=0; $t<@to_process; ++$t )
    {
    $vals{$to_process[$t]} =~ s/!//g;
    }

  $positions = join(",", @positions);
  $positions =~ s/!//g;
  $positions =~ s/,{2,}/,/g;
  $positions =~ s/^,//g;
  $vals{'POSITIONS'} = $positions;

  return(%vals);
  }
return 1;  