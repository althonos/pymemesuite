#!/usr/bin/perl
use strict;
use warnings;
# Replace the key words "@maxtime@" and "@maxtime_short@" 
# with the times, converted from seconds to hh:mm:ss format.
my $maxtime = $ARGV[0]; 	# maxtime
my $maxtime_short =$ARGV[1]; 	# maxtime_short
my ($i, $h, $m, $s);
while (<STDIN>) {
  for ($i=0; $i<2; $i++) {
    my $type = ($i==0) ? "maxtime" : "maxtime_short";
    my $t = ($i==0) ? $maxtime : $maxtime_short;
    $h=int($t/3600); 
    $m=int(($t-$h*3600)/60); 
    $s=$t-$h*3600-$m*60;
    $m = "0$m" if ($m < 10);
    $s = "0$s" if ($s < 10);
    s/\@$type\@/$h:$m:$s/g; 
  }
  print;
}
