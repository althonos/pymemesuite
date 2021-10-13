package MotifInStremeXML;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

use StremeSAX;

sub new {
  my $classname = shift;
  my $self = {};
  bless($self, $classname);
  $self->_init(@_);
  return $self;
}

sub parse_more {
  my $self = shift;
  my ($buffer) = @_;
  $self->{parser}->parse_more($buffer);
}

sub parse_done {
  my $self = shift;
  $self->{parser}->parse_done();
}

sub has_errors {
  my $self = shift;
  return ($self->{parser}->has_errors() or $self->{alph}->has_errors());
}

sub get_errors {
  my $self = shift;
  return ($self->{parser}->get_errors(), $self->{alph}->get_errors());
}

sub get_motifs {
  my $self = shift;
  my @motifs = @{$self->{queue}};
  $self->{queue} = [];
  return @motifs;
}

sub has_motif {
  my $self = shift;
  return scalar(@{$self->{queue}}) > 0;
}

sub next_motif {
  my $self = shift;
  return pop(@{$self->{queue}});
}

sub _init {
  my $self = shift;

  my $sax = new StremeSAX($self, 
    'handle_alphabet' => \&_handle_alphabet,
    'handle_strands' => \&_handle_strands,
    'handle_background' => \&_handle_background,
    'handle_test_positives' => \&_handle_test_positives,
    'start_motif' => \&_start_motif,
    'end_motif' => \&_end_motif,
    'handle_pos' => \&_handle_pos
  );
  $self->{parser} = $sax;
  $self->{queue} = [];
  $self->{motif} = undef;
  $self->{alph} = undef;
  $self->{strands} = undef;
  $self->{bg} = {};
  $self->{test_positives} = undef;

}

sub _handle_alphabet {
  my ($self, $alphabet) = @_;
  $self->{alph} = $alphabet;
  $self->{bg}->{alph} = $self->{alph};
}

sub _handle_strands {
  my $self = shift;
  my ($strands) = @_;
  $self->{strands} = ($strands eq 'both' ? 2 : ($strands eq 'forward' ? 1 : 0));
}

sub _handle_background {
  my $self = shift;
  my (@probs) = @_;
  my $bg = $self->{bg};
  my $alph = $self->{alph};
  for (my $i = 0; $i < $alph->size_core(); $i++) {
    my $sym = $alph->char($i);
    $bg->{$sym} = $probs[$i];
  }
}

sub _handle_test_positives{
  my ($self, $test_positives) = @_;
  $self->{test_positives} = $test_positives;
}

sub _start_motif {
  my $self = shift;
  my ($id, $alt, $width, $initial_width, $seed, $score_threshold, 
    $train_pos_count, $train_neg_count, $train_log_pvalue, $train_pvalue, $train_dtc, $train_bernoulli,
    $test_pos_count, $test_neg_count, $test_log_pvalue, $test_pvalue, $test_log_evalue, $test_evalue,
    $test_dtc, $test_bernoulli, $elapsed_time, $total_sites, $site_distr) = @_;

  my %pspm = ();
  my $alph = $self->{alph};
  for (my $i = 0; $i < $alph->size_core(); $i++) {
    my $sym = $alph->char($i);
    $pspm{$sym} = [];
  }

  my $evalue_accurate = ($self->{test_positives} > 0);
  my $motif = {
    bg => $self->{bg},
    strands => $self->{strands},
    id => $id, 
    alt => $alt,
    width => $width, 
    sites => $train_pos_count + $test_pos_count,
    pseudo => 0,
    log_evalue_accurate => $evalue_accurate,
    log_evalue => $evalue_accurate ? $test_log_evalue : $train_log_pvalue,
    evalue => $evalue_accurate ? $test_evalue : $train_pvalue,
    log_pvalue => $evalue_accurate ? $test_log_pvalue : $train_log_pvalue,
    pspm => \%pspm
  };

  $self->{motif} = $motif;
}

sub _end_motif {
  my $self = shift;
  push(@{$self->{queue}}, $self->{motif});
  $self->{motif} = undef;
}

sub _handle_pos {
  my $self = shift;
  my (@probs) = @_;
  my $alph = $self->{alph};
  my $pspm = $self->{motif}->{pspm};
  for (my $i = 0; $i < $alph->size_core(); $i++) {
    my $sym = $alph->char($i);
    push(@{$pspm->{$sym}}, $probs[$i]);
  }
}

1;
