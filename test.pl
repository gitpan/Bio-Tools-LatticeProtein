use Test;
BEGIN { plan tests => 5 }
ok(1);

use Data::Dumper;
use List::Util qw(min);
use Bio::Tools::LatticeProtein;

$lp = Bio::Tools::LatticeProtein->new(
				      sequence => 'perlisnifty',
				      cycle => 15,
				      mutdist => {
					  end_rotation => 0.1,
					  kink_jump    => 0.4,
					  crankshaft   => 0.4,
					  slithering   => 0.1,
				      },
				      );
$lp->run;

ok($lp->{sequence}, 'PERLISNIFTY');
ok({qw,-2 1 -3 1,}->{min(map{$lp->{pool}->[$_]->{energy}}0..$lp->{popsize})}, 1);
ok($lp->{mutdist}->{kink_jump}, 0.5);

ok($lp->result);

