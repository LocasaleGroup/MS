#!/usr/bin/perl
#  Copyright Hunter Moseley, 2009. All rights reserved.
#  Written by Hunter Moseley 7/6/2009

my $help = <<HELP;
carbon_na_stripping.pl - strip carbon natural abundance contribution from given mass isotopomer data.

 usage: carbon_na_stripping.pl carbon_count [-excel] < data_sets

 Parameters:
       carbon_count - number of carbons in the detected molecular entity.
       -excel - print out a comma delimited output with no extra columns besides data and stripped.
       data_set - space or comma delimited list of mass isotopomers for each 
HELP

use FindBin;
use lib $FindBin::Bin;

use strict;

if (! @ARGV || $ARGV[0] =~ /^-h/)
  {
  print STDERR $help;
  exit;
  }

my $max_carbons = shift @ARGV;
die "Error:  zero or negative carbon count \"$max_carbons\" given." if ($max_carbons <= 0);

my $excel = 0;
if (@ARGV && $ARGV[0] =~ /^-e/)
  { $excel = 1; }

# calculate carbon NA intermediate values needed to strip NA contributions out.
my @pure_na = (1);
for(my $x=1; $x <= $max_carbons; $x++)
  { $pure_na[$x] = $pure_na[$x-1] * 0.01109; }

my @pure_reverse_na = (1);
for(my $x=1; $x <= $max_carbons; $x++)
  { $pure_reverse_na[$x] = $pure_reverse_na[$x-1] * (1 - 0.01109); }


my @binomial_na;
my @binomial_na_sum;
for(my $n=0; $n < $max_carbons; $n++)
  {
  for(my $k=$max_carbons; $k > $n; $k--)
    { 
    $binomial_na[$n][$k] = &binomial($max_carbons-$n,$k-$n) * $pure_na[$k-$n] * $pure_reverse_na[$max_carbons-$k]; 
    $binomial_na_sum[$n] += $binomial_na[$n][$k];
    }
  }

my @excel_columns;
# process each dataset
foreach my $line (<STDIN>)
  {
  chomp $line;
  if ($line =~ /^\#/) # skip commented lines
    { 
    print "$line\n"; 
    next; 
    }
  next if ($line =~ /^\s*$/); # skip empty lines
  $line =~ tr/[,]/ /;
  my $data = [ split(/\s+/,$line) ];
  my $data_size = scalar(@$data);
  my ($converted,$calc, $num_iter, $negative) = &remove_NA($data);
  my $full_diff = 0; map { $full_diff += abs($$data[$_] - $$calc[$_]); } (0..$#$calc);
  my $obs_diff = 0; map { $obs_diff += abs($$data[$_] - $$calc[$_]); } (0..$data_size);
  my @missing = map { $$calc[$_]; } (grep { $$data[$_] == 0; } (0..$#$calc));
  my $missing_sum = &sum(@missing);
  my $largest_real = (sort { $b <=> $a } (@$data))[0];
  my $largest_missing = 0;
  if (@missing)
    { $largest_missing = (sort { $b <=> $a; } (@missing))[0]; }
  my $real_data_sum = &sum(@$data);
  my $data_sum = &sum(map { ($$data[$_] <= 0) ? $$calc[$_] : $$data[$_]; } (0..$#$calc));
  my $converted_sum = &sum(@$converted);
  my $renormalized = [ map { $_ / $converted_sum * $data_sum; } (@$converted) ];
  if (! $excel)
    {
    print "Data: @$data\n";
    print "Calculated: @$calc\n";
    print "Num iterations: \t$num_iter\n";
    print "Full Difference: \t$full_diff / $data_sum = ",$full_diff / $data_sum,"\n";
    print "Observable Difference: \t$obs_diff / $real_data_sum = ",$obs_diff / $real_data_sum,"\n";
    print "Fraction Negative: \t$negative / $real_data_sum = ", $negative / $real_data_sum,"\n";
    print "Fraction Missing: \t$missing_sum / $real_data_sum = ", $missing_sum / $real_data_sum, "\n";
    print "Largest Missing: \t$largest_missing / $largest_real = ", $largest_missing / $largest_real, "\n";
    print "NA Stripped: @$converted\n";
    print "Renormalized NA Stripped: @$renormalized\n";
    }
  else
    { push @excel_columns,$data, $calc, $renormalized; }
  }

if ($excel)
  {
  print "\nIsotopologue, ","Data, Calc, Stripped, "x(@excel_columns/3),"\n";
  my $rows = scalar(@{$excel_columns[$#excel_columns]});
  for(my $x=0; $x < $rows; $x++)
    { print "$x, ",join(", ", map { $$_[$x] * 1; } (@excel_columns)),"\n"; }
  }

exit;


sub sum
  {
  my $sum = 0;
  foreach my $value (@_)
    { $sum += $value; }

  return $sum;
  }

sub binomial
  {
  my $n = shift @_;
  my $k = shift @_;
  my $diff = $n - $k;

  my $value = 1;

  if ($k > $diff)
    {
    while($n > $k)
      { $value *= $n; $n--; }
    
    while($diff > 0)
      { $value /= $diff; $diff--; }
    }
  else
    {
    while($n > $diff)
      { $value *= $n; $n--; }
    
    while($k > 0)
      { $value /= $k; $k--; }
    }

  return $value;
  }

sub remove_NA
  {
  my $data = shift @_;
  my $num_iter = 1;

  my $target_sum = &sum(@$data);
  my ($subtracted, $negative) = &subtract_NA($data);
  my $subtracted_sum = &sum(@$subtracted);
  my $added = &add_NA([ map { $_ * $target_sum / $subtracted_sum; } (@$subtracted) ]);
  $target_sum = &sum(map { ($$data[$_] <= 0) ? $$added[$_] : $$data[$_]; } (0..$#$added));
  my $added_sum = &sum(@$added);
  
  my $diff = 0; map { $diff += abs($$data[$_] - $$added[$_]); } (0..$#$added);
  my $last_diff = $diff * 2; 
  
  while($diff < $last_diff)
    {
    $num_iter++;
    $last_diff = $diff;
    my $data2 = [ map { ($$data[$_] <= 0) ? $$added[$_] : $$data[$_]; } (0..$#$added) ];
    my $data2_sum = &sum(@$data2);
    ($subtracted, $negative) = &subtract_NA([ map { $_ * $target_sum / $data2_sum; } (@$data2) ]);
    $subtracted_sum = &sum(@$subtracted);
    $added = &add_NA([ map { $_ * $target_sum / $subtracted_sum; } (@$subtracted) ]);
    $added_sum = &sum(@$added);
    $target_sum = &sum(map { ($$data[$_] <= 0) ? $$added[$_] : $$data[$_]; } (0..$#$added));
    $diff = 0; map { $diff += abs($$data[$_] - $$added[$_]); } (0..$#$added);
    }
  
  return ($subtracted, $added,$num_iter, $negative);
  }

sub subtract_NA
  {
  my $orig = shift @_;
  my $calc = [ @$orig ];
  my $negative = 0;
  
  for(my $y=0; $y <= $max_carbons; $y++)
    {
    my $x=0;
    while($x < $y)
      { 
      $$calc[$y] -= $$calc[$x] * $binomial_na[$x][$y]; 
      $x++;
      }    

    $$calc[$y] /= (1 - $binomial_na_sum[$y]);
    if ($$calc[$y] < 0) # must do this to improve accuracy.
      {
      $negative += $$calc[$y];
      $$calc[$y] = 0; 
      }
    }

  return ($calc, $negative);
  }


sub add_NA
  {
  my $orig = shift @_;
  my $calc = [ @$orig ];

  for(my $y=$max_carbons; $y >= 0; $y--)
    {
    $$calc[$y] *= (1 - $binomial_na_sum[$y]);
    my $x=0;
    while($x < $y)
      { 
      $$calc[$y] += $$calc[$x] * $binomial_na[$x][$y]; 
      $x++;
      } 
    }

  return $calc;
  }