package StremeSAX;

use base 'CheckingSAX';

use strict;
use warnings;

use CheckingSAX;

my $num_re = qr/^((?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|inf)$/;
my $num_trim_re = qr/^\s*((?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|inf)\s*$/;
my $int_re = qr/^(\d+)$/;
my $pos_int_trim_re = qr/^\s*(\d+)\s*$/;
my $int_trim_re = qr/^\s*([-+]?\d+)\s*$/;
my $float_re = qr/^(\d+(?:\.\d+)?)$/;
my $float_trim_re = qr/^\s*(\d+(?:\.\d+)?)\s*$/;
my $trim_re = qr/^\s*(.*?)\s*$/;
my $text_re = qr/^(.*\S.*)$/;
my $nospace_re = qr/^(\S+)$/;
my $allowed_sym_re = qr/[!-~]/; # Allow all visible ASCII
my $sym_re = qr/^($allowed_sym_re)$/;
my $syms_re = qr/^($allowed_sym_re+)$/;
my $sym_id_re = qr/^([a-zA-Z]|n[0-9]|x[a-zA-Z0-9]{2})$/;
my $num_list_re = qr/^\s*((\d+\s+)*(\d+))\s*$/;

my $ST_IN_MOTIF = [
  {ELE => 'pos', RPT => $RPT_ALO, ATRS => [\&_get_alphabet_atrs]}
];

my $ST_IN_MOTIFS = [
  {ELE => 'motif', RPT => $RPT_ANY, TO => $ST_IN_MOTIF, ATRS => [
      {ATR => 'id', VAL => $nospace_re},
      {ATR => 'alt', VAL => $nospace_re},
      {ATR => 'width', VAL => $int_re},
      {ATR => 'initial_width', VAL => $int_re},
      {ATR => 'seed', VAL => $nospace_re},
      {ATR => 'score_threshold', VAL => $num_re},
      {ATR => 'train_pos_count', VAL => $int_re},
      {ATR => 'train_neg_count', VAL => $int_re},
      {ATR => 'train_log_pvalue', VAL => $num_re},
      {ATR => 'train_pvalue', VAL => $num_re},
      {ATR => 'train_dtc', VAL => $num_re},
      {ATR => 'train_bernoulli', VAL => $num_re},
      {ATR => 'test_pos_count', VAL => $int_re},
      {ATR => 'test_neg_count', VAL => $int_re},
      {ATR => 'test_log_pvalue', VAL => $num_re},
      {ATR => 'test_pvalue', VAL => $num_re},
      {ATR => 'test_log_evalue', VAL => $num_re},
      {ATR => 'test_evalue', VAL => $num_re},
      {ATR => 'test_dtc', VAL => $num_re},
      {ATR => 'test_bernoulli', VAL => $num_re},
      {ATR => 'elapsed_time', VAL => $num_re},
      {ATR => 'total_sites', OPT => $int_re},
      {ATR => 'site_distr', OPT => $num_list_re},
      {ATR => 'max_sites', OPT => $int_re},
      {ATR => 'site_hist', OPT => $num_list_re}
    ]
  }
];

my $ST_IN_MODEL = [
  {ELE => 'command_line', VAL => $text_re},
  {ELE => 'train_positives', ATRS => [
      {ATR => 'count', VAL => $int_re},
      {ATR => 'positions', VAL => $int_re},
      {ATR => 'maxlen', OPT => $int_re},
      {ATR => 'file'}
    ]},
  {ELE => 'train_negatives', ATRS => [
      {ATR => 'count', VAL => $int_re},
      {ATR => 'positions', VAL => $int_re},
      {ATR => 'from', VAL => qr/^(file|shuffled|none)$/},
      {ATR => 'file', OPT => $text_re}
    ], WATCH => \&_watch_ele_negatives},
  {ELE => 'test_positives', ATRS => [
      {ATR => 'count', VAL => $int_re},
      {ATR => 'positions', VAL => $int_re}
    ]},
  {ELE => 'test_negatives', ATRS => [
      {ATR => 'count', VAL => $int_re},
      {ATR => 'positions', VAL => $int_re}
    ]},
  &_struc_alphabet(ATRS_AS => 'alphabet_atrs', ALPH_AS => 'alph', ID2SYM_AS => 'id_to_sym', SYM2ID_AS => 'sym_to_id'),
  {ELE => 'strands', VAL => qr/^(both|given|none)$/},
  {ELE => 'sequence_db', ATRS => [
      \&_get_alphabet_atrs
    ]},
  &_struc_alph_freqs('background_frequencies', 
    [ {ATR => 'source', VAL => $text_re},
      {ATR => 'order', VAL => $pos_int_trim_re}
    ]),
  {ELE => 'stop', ATRS => [
      #{ATR => 'pvt', OPT => \&_value_pvalue}, 
      {ATR => 'thresh_type', OPT => $trim_re}, 
      {ATR => 'thresh', OPT => \&_value_evalue}, 
      {ATR => 'nmotifs', OPT => $int_re}, 
      {ATR => 'time', OPT => $num_trim_re}
  ]},
  {ELE => 'objfun', VAL => $trim_re},
  {ELE => 'test', VAL => $trim_re},
  {ELE => 'minw', VAL => $pos_int_trim_re},
  {ELE => 'maxw', VAL => $pos_int_trim_re},
  {ELE => 'kmer', VAL => $pos_int_trim_re},
  {ELE => 'hofract', VAL => $num_trim_re},
  {ELE => 'neval', VAL => $pos_int_trim_re},
  {ELE => 'nref', VAL => $pos_int_trim_re},
  {ELE => 'niter', VAL => $pos_int_trim_re},
  {ELE => 'patience', VAL => $pos_int_trim_re},
  {ELE => 'seed', VAL => $pos_int_trim_re},
  {ELE => 'useer', VAL => \&_value_yesno},
  {ELE => 'minscore', VAL => $int_trim_re},
  {ELE => 'ignore_depth', VAL => $int_trim_re},
  {ELE => 'nsubsets', VAL => $int_trim_re},
  {ELE => 'min_pal_ratio', VAL => $num_trim_re},
  {ELE => 'max_pal_ed', VAL => $num_trim_re},
  {ELE => 'cand', VAL => \&_value_yesno},
  {ELE => 'experimental', VAL => \&_value_yesno},
  {ELE => 'totallength', VAL => $pos_int_trim_re},
  {ELE => 'align', VAL => $trim_re},
  {ELE => 'host', VAL => $trim_re},
  {ELE => 'description', RPT => $RPT_OPT, VAL => qr/^([\S\s]*)$/}
];
my $ST_IN_STREME = [
  {ELE => 'model', TO => $ST_IN_MODEL}, 
  {ELE => 'motifs', TO => $ST_IN_MOTIFS},
  {ELE => 'reason_for_stopping', VAL => $trim_re},
  {ELE => 'run_time', ATRS => [
      {ATR => 'cpu', VAL => $trim_re} ]}
];
my $ST_START = [
  {ELE => 'STREME', NAME => 'streme', RPT => $RPT_ONE, TO => $ST_IN_STREME, ATRS => [
      {ATR => 'version', VAL => qr/^(\d+)\.(\d+)\.(\d+)$/, EXPAND => 3},
      {ATR => 'release', VAL => $text_re}
    ]}];

#
# Return the definition of a STREME XML file.
sub _xml_def {
  return $ST_START;
}

sub _init {
  my $self = shift;
  $self->{alphabet_atrs} = [];
  $self->SUPER::_init(@_);
}

sub _get_alphabet_atrs {
  return @{$_[0]->{alphabet_atrs}};
}

sub _watch_ele_negatives {
  my ($self, $attrs, $name, $count, $from, $file, $last_mod) = @_;
  if ($from eq 'file') {
    $self->_error('negatives/@file missing') unless defined $file;
  }
}

sub _struc_alph_array {
  my ($callback, %options) = @_;
  my %values;
  my $fn_init = sub {
    %values = ();
  };
  my $fn_letter_value = sub {
    my ($self, $id, $value) = @_;
    if (defined($values{$id})) {
      $self->_error("Value repeated for letter ID \"$id\".");
    }
    $values{$id} = $value;
  };
  my $fn_result = sub {
    my ($self) = @_;
    my @alph_array = ();
    for (my $i = 0; $i < $self->{alph}->size_core(); $i++) {
      my $sym = $self->{alph}->char($i);
      my $id = $self->{sym_to_id}->{$sym};
      my $value = $values{$id};
      unless (defined($value)) {
        $self->_error("Value missing for letter ID \"$id\".");
        $value = 0;
      }
      $alph_array[$i] = $value;
    }
    $callback->(@alph_array);
  };
  return (
    {
      ELE => 'alphabet_array',
      RPT => $RPT_ALO,
      SILENT => 1,
      START => $fn_init,
      END => $fn_result,
      TO => [
        {
          ELE => 'value',
          RPT => $RPT_ALO,
          SILENT => 1,
          END => $fn_letter_value,
          ATRS => [
            {ATR => 'letter_id', VAL => $sym_id_re}
          ],
          VAL => $num_trim_re
        }
      ]
    }
  );
}
sub _struc_alph_freqs {
  my ($element, $atrs) = @_;
  $atrs = [] unless defined $atrs;
  my @freqs;
  my $fn_array = sub {
    @freqs = @_;
  };
  my $fn_result = sub {
    my ($self, @atrs) = @_;
    return (@atrs, \@freqs);
  };
  return (
    {
      ELE => $element,
      HANDLE => $fn_result,
      TO => [&_struc_alph_array($fn_array)],
      ATRS => $atrs
    }
  );
}
