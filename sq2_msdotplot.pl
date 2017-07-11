#!/usr/bin/perl
#
# A perl script to generate variation dotplots as described in;
#
# http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006862
#
# Please see README and LICENSE for details of use
#
#=========================================================================================

use Text::Diff;

$filein = @ARGV[0];
$gene = $filein;
$filein =~ tr/A-Z/a-z/;
chomp($filein);
$fileout = $filein;

$scale_level = @ARGV[1];
chomp($scale_level);

$colour = @ARGV[2];
chomp($colour);


print scalar(localtime(time))," Start Time\n";

$filein = join("", 'libs/', $filein, '.exon23.lib');

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
  ($allele,$seq) = split(/\t/, $_);

  print "Using $allele as reference sequence for INDEL removal\n";
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
  ($allele,$seq) = split(/\t/, $_);

  $seq{$allele} = $seq;
  }
close(LIB);

@allele = keys(%seq);
for ( $a=0; $a<@allele; ++$a )
  {
  $allele[$a] =~ s/\_/\*/;  
  }
@allele = sort_new_nomenc::sorting_nomenclature(@allele);
for ( $a=0; $a<@allele; ++$a )
  {
  $allele[$a] =~ s/\*/\_/;  
  }

# Remove allele with indels

for ( $bp=0; $bp<length($refseq); ++$bp )
  {
  $refbase = substr($refseq, $bp, 1);  
  if ( $refbase eq '.' )
    {
    substr($refseq, $bp, 1, '!');
    for ( $a=0; $a<@allele; ++$a )
      {
      $seqbase =  substr($seq{$allele[$a]}, $bp, 1); 
      if ( $seqbase eq '.' )  
        {
        substr($seq{$allele[$a]}, $bp, 1, '!');
        }
      elsif( $seqbase =~ /[ACTGX]/ )  
        {
        substr($seq{$allele[$a]}, $bp, 1, 'I');
        }
      }
    }
  }

for ( $a=0; $a<@allele; ++$a )
  {
  $seq{$allele[$a]} =~ s/I{1,}[ACGTX]/I/g;  
  $seq{$allele[$a]} =~ s/\./D/g;  
  $seq{$allele[$a]} =~ s/!//g;  
  push(@include, $allele[$a]);
  }
  
@allele = @include;

# Substitute unsequenced regions

&define_substitutions;

for ( $a=0; $a<@allele; ++$a )
  {
  if ( $substitute{$allele[$a]} =~ /:/ && $seq{$allele[$a]} =~ /\*/ )
    {
    print "Sequence of $allele[$a] contains unsequenced region. Substituting $substitute{$allele[$a]} for $allele[$a] where needed\n";
    print EXCLUDE "Substituting $substitute{$allele[$a]} for $allele[$a] where needed\n";
    $sub_allele = $substitute{$allele[$a]};
    push(@substitutions, $allele[$a]);
    for ( $bp=0; $bp<length($seq{$allele[$a]}); ++$bp )
      {
      $seqbase = substr($seq{$allele[$a]}, $bp, 1);  
      if ( $seqbase eq '*' )
        {
        $newbase =  substr($seq{$sub_allele}, $bp, 1); 
        substr($seq{$allele[$a]}, $bp, 1, $newbase);
        }
      }
    }
  if ( $seq{$allele[$a]} =~ /\*/ )
    {
    print "WARNING! Unsequenced region remains in $allele[$a]\n";
    }  
  }
  

# Remove monomorphic positions

for ( $bp=0; $bp<length($seq{$allele[0]}); ++$bp )
  {
  $refbase = substr($seq{$allele[0]}, $bp, 1);  
  $mono = 1;
  for ( $a=1; $a<@allele; ++$a )
    {
    $seqbase = substr($seq{$allele[$a]}, $bp, 1);  
    if ( $seqbase ne $refbase )
      {
      $pos = $bp+101;  
      $mono = 0;
      last;
      }
    }
  if ( $mono eq 1 )
    {
    for ( $a=0; $a<@allele; ++$a )
      {
      substr($seq{$allele[$a]}, $bp, 1, '!');
      }
    }
  }
  
for ( $a=0; $a<@allele; ++$a )
  {
  $seq{$allele[$a]} =~ s/!//g;
  }


@allele = @include;
@allele = sort_new_nomenc::sorting_nomenclature(@allele);
$no_alleles = @allele;

print "Running on $no_alleles alleles\n";

$max_num_mismatch = 0;

use GD;
$font = "./Vera.ttf";

$image = new GD::Image(250+$no_alleles,250+$no_alleles);
&define_colours;

$image->interlaced('true');
$image->rectangle(0,0,150+$no_alleles,175+$no_alleles,$white);

# Print Dot Plot

$adjust = 20*($no_alleles/10000);
$adjust = sprintf "%.0f", $adjust;

print "Axis size of $adjust\n";

$image->filledRectangle(99-$adjust,99,99,99+$no_alleles,$black);
$image->filledRectangle(99-$adjust,100+$no_alleles,100+$no_alleles,100+$adjust+$no_alleles,$black);

$x_group = '00';
$y_group = '00';

for ( $a=0; $a<@allele; ++$a )
  {
  $x_new_group = substr($allele[$a], 2, 2);
  if ( $x_group ne $x_new_group )
    {
    $image->filledRectangle(75,98+$a,100,102+$a,$black);
    $x_group = $x_new_group;   
    }
  }

$prev_b = 0;
for ( $b=0; $b<@allele; ++$b )
  {
  $y_new_group = substr($allele[$b], 2, 2);
  if ( $y_group ne $y_new_group )
    {
    $image->filledRectangle(98+$b,100+$no_alleles,102+$b,125+$no_alleles,$black);
    $y_group = $y_new_group;
      
    $b_size = $b-$prev_b;
    $b_size = 'N = ' . $b_size;
    $prev_b = $b;  
    }  
  }
  
  
for ( $a=0; $a<@allele; ++$a )
  {
  print "Parsing $allele[$a]\n";
  for ( $b=$a; $b<@allele; ++$b )
    {
    $num_mismatch = hd($seq{$allele[$a]}, $seq{$allele[$b]});
    if ( $scale_level < 1 )
      {
      $scale_level = 1;  
      }
    
    $num_mismatch = $num_mismatch/$scale_level;
    
    if ( $num_mismatch < 21 )
      {
      $image->setPixel(100+$a,100+$b,$scale[$num_mismatch]);
      }
    elsif ( $num_mismatch > 20 )
      {
      $image->setPixel(100+$a,100+$b,$scale[21]);  
      }
    }
  }

# print Scale

for ( $s=0; $s<@scale; ++$s )
  {
  if ( $s eq 21 )
    {
    $image->filledRectangle(100+$s*$adjust*2+$s*$adjust*8,150+$no_alleles,100+$s*$adjust*2+($s+1)*$adjust*8,150+$no_alleles+$adjust*8,$black);
    $image->filledRectangle(102+$s*$adjust*2+$s*$adjust*8,152+$no_alleles,98+$s*$adjust*2+($s+1)*$adjust*8,148+$no_alleles+$adjust*8,$white);
    }
  else
    {
    $image->filledRectangle(100+$s*$adjust*2+$s*$adjust*8,150+$no_alleles,100+$s*$adjust*2+($s+1)*$adjust*8,150+$no_alleles+$adjust*8,$scale[$s]);
    
    print "    $image->filledRectangle(100+$s*$adjust*2,150+$no_alleles,100+$s*$adjust*2+$adjust*8,150+$no_alleles+$adjust*8,$scale[$s]);\n";
    }  
  }
  
$png_data = $image->png;
system("rm temp.png");
open (DISPLAY,">temp.png") || die;
binmode DISPLAY;
print DISPLAY $png_data;
close DISPLAY;

print "Based on $no_alleles alleles\n";  
exit;

sub hd
  {
  return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
  }

sub define_substitutions  
  {
  $substitute{'A_31:57'} = 'A_31:01:02';
  
  $substitute{'B_07:08'} = 'B_07:02:01';
  $substitute{'B_08:03'} = 'B_08:01:01';
  $substitute{'B_45:05'} = 'B_45:01:01';
  $substitute{'B_51:10'} = 'B_51:01:01';
  $substitute{'B_35:167'} = 'B_35:168';

  $substitute{'C_07:65'} = 'C_07:60';
  $substitute{'C_12:60'} = 'C_12:02:01';
  }
  
sub define_colours
  {
  $white     = $image->colorAllocate(255,255,255);
  $black     = $image->colorAllocate(0,0,0);       
  
  if ( $colour eq 1 )
    {
    $scale[0]  = $image->colorAllocate(213,0,2);
  
    $scale[1]  = $image->colorAllocate(212,34,0);
    $scale[2]  = $image->colorAllocate(211,70,0);
    $scale[3]  = $image->colorAllocate(210,106,0);
    $scale[4]  = $image->colorAllocate(209,142,0);
    $scale[5]  = $image->colorAllocate(208,177,0);
    $scale[6]  = $image->colorAllocate(202,207,0);
    $scale[7]  = $image->colorAllocate(166,206,0);
    $scale[8]  = $image->colorAllocate(130,205,0);
    $scale[9]  = $image->colorAllocate(94,204,0);
    $scale[10] = $image->colorAllocate(58,203,0);
  
    $scale[11] = $image->colorAllocate(23,202,0);
    $scale[12] = $image->colorAllocate(0,202,11);
    $scale[13] = $image->colorAllocate(0,201,45);
    $scale[14] = $image->colorAllocate(0,200,79);
    $scale[15] = $image->colorAllocate(0,199,113);
    $scale[16] = $image->colorAllocate(0,198,147);
    $scale[17] = $image->colorAllocate(0,197,180);
    $scale[18] = $image->colorAllocate(0,179,196);
    $scale[19] = $image->colorAllocate(0,144,195);
    $scale[20] = $image->colorAllocate(0,110,194);
  
    $scale[21] = $image->colorAllocate(255,255,255);
    }
  else # Greyscale 
    {
    $scale[0]  = $image->colorAllocate(0,0,0);
   
    for ( $rgb=1; $rgb<6; ++$rgb)
      {
      $scale[$rgb]  = $image->colorAllocate($rgb*20,$rgb*20,$rgb*20);      
      }
    for ( $rgb=1; $rgb<16; ++$rgb)
      {
      $scale[$rgb+5]  = $image->colorAllocate($rgb*10+100,$rgb*10+100,$rgb*10+100);      
      }
   
    $scale[21] = $image->colorAllocate(255,255,255);
    }
  }
  
