package Bio::Tools::LatticeProtein::Vars;

use strict;
our $VERSION = '0.01';

use Exporter;
our @EXPORT = qw(%hp_type @dir);
our @ISA = qw(Exporter);

our %hp_type = qw/
  A 1
  C 2
  D 2
  E 2
  F 1
  G 2
  H 2
  I 1
  K 2
  L 1
  M 1
  N 2
  P 1
  Q 2
  R 2
  S 2
  T 2
  V 1
  W 1
  Y 2
    /;

our @dir = (
            [ qw/0 0 1/ ],
            [ qw/0 1 0/ ],
            [ qw/1 0 0/ ],
            [ qw/0 0 -1/ ],
            [ qw/0 -1 0/ ],
            [ qw/-1 0 0/ ],
            );
1;
