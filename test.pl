use Test;
BEGIN { plan tests => 6 }
ok(1);

use Data::Dumper;
use List::Util qw(min);
use Bio::Tools::LatticeProtein;

$lp = Bio::Tools::LatticeProtein->new(
				      sequence => 'perlisnifty',
				      cycle => 100,
				      mutdist => {
					  end_rotation => 0.2,
					  kink_jump    => 0.3,
					  crankshaft   => 0.3,
					  slithering   => 0.2,
				      },
				      );

$lp->create_pool;

ok($lp->is_continuous(0), 1);

$lp->run;

ok($lp->{sequence}, 'PERLISNIFTY');

#print min(map{$lp->{pool}->[$_]->{energy}}0..$lp->{popsize});

ok({map{-$_=>1} 2..7}->{min(map{$lp->{pool}->[$_]->{energy}}0..$lp->{popsize})}, 1);
ok($lp->{mutdist}->{kink_jump}, 0.5);

ok($lp->result);

