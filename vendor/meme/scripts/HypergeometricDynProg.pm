# AUTHOR: Timothy L. Bailey
# CREATE DATE: 7-5-2005
# MODIFIER: Philip Machanick
# MODIFIED DATE: 9-08-2009
# changed to use dynamic programming version of combinations
package HypergeometricDynProg;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(hypergeometric_dp);

use CombDynProg qw(log_comb_dp);

################################################################################
#        hypergeometric
#
#        
#        The probability of finding k red balls among n draws without
#        replacement from an urn with K red balls and N-K black ones.
#        
#
################################################################################
sub hypergeometric_dp{
  my(
    $k,                                # number of red balls drawn
    $n,                                # number of draws
    $K,                                # number of red balls in urn
    $N                                # number of balls in urn
  ) = @_;
  my $p;
  
  if ($k > $K) {
    $p = 0;
  } else {
    $p = exp(log_comb_dp($K, $k) + log_comb_dp($N-$K, $n-$k) - log_comb_dp($N, $n));
  }
  return($p);
} # hypergeometric

