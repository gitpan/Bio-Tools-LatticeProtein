package Bio::Tools::LatticeProtein;

use strict;
our $VERSION = '0.01';
use Bio::Tools::LatticeProtein::Vars;
use List::Util qw(shuffle);

######################################################################
# Initializing protein pool
######################################################################

sub create_fold {
    my $pkg = shift;
    my $foldidx = shift;
    my @chain;
    my %idx;

    my $cnt = 0;
    foreach my $i (0..$pkg->{length}-1){
	$idx{join(q/,/, $cnt, 0, 0)} = $pkg->{typearr}->[$i];
	push @chain, [ $cnt, 0, 0 ];
	$cnt++;
    }

    $pkg->{pool}->[$foldidx] = { chain => \@chain, idx => \%idx, energy => 0 };
}

sub create_pool {
    my $pkg = shift;
    for( 0..$pkg->{popsize}-1){
	$pkg->create_fold($_);
    }
}


######################################################################
# Constructor
######################################################################

sub new {
    my $pkg = shift;
    my %arg = @_;

    die "Survival should not exceeds Population size\n" if $arg{survival} > $arg{popsize};
    die "The sum of mutation probability distribution must be equal to ONE\n" if
	defined $arg{mutdist} && 
	$arg{mutdist}->{end_rotation} +
	    $arg{mutdist}->{kink_jump} +
		$arg{mutdist}->{crankshaft} +
		    $arg{mutdist}->{slithering} != 1;
    my $def_mutdist = {
	end_rotation => 0.25,
	kink_jump    => 0.5,
	crankshaft   => 0.75,
	slithering   => 1,
    };

    if(defined $arg{mutdist}){
	$arg{mutdist}->{kink_jump}  += $arg{mutdist}->{end_rotation};
	$arg{mutdist}->{crankshaft} += $arg{mutdist}->{kink_jump};
	$arg{mutdist}->{slithering} = 1;
    }


    bless {
        sequence   => uc($arg{sequence} || die "Please input a sequence\n"),
	length     => length($arg{sequence}),
	typearr    => [ map { $hp_type{$_} } split( //, uc $arg{sequence}) ],

        cycle      => ($arg{cycle}      || 15),
	popsize    => ($arg{popsize}    || 20),
	survival   => ($arg{survival}   || 10),
	mutdist    => (defined $arg{mutdist} ? $arg{mutdist} : $def_mutdist),

	pool       => undef,
    }, $pkg;
}


######################################################################
# Utils
######################################################################
sub one_step    { [ $_[0]->[0]+$_[1]->[0], $_[0]->[1]+$_[1]->[1], $_[0]->[2]+$_[1]->[2] ] }
sub is_pos_free { !$_[0]->{idx}->{join q/,/, @{$_[1]} } }
sub pos_equal {$_[0]->[0]==$_[1]->[0] && $_[0]->[1]==$_[1]->[1] && $_[0]->[2]==$_[1]->[2]}
sub possub { [ map{$_[0]->[$_]-$_[1]->[$_]} 0..2 ] }
sub posadd { [ map{$_[0]->[$_]+$_[1]->[$_]} 0..2 ] }
sub select_freepos {
    my $fold = shift;
    my $pos = shift;
    foreach my $dir (shuffle @dir){
	my $new_pos = one_step($pos, $dir);
	return $new_pos if is_pos_free($fold, $new_pos);
    }
    undef;
}

sub is_dist_one {
    my @dv = map { abs($_[0]->[$_] - $_[1]->[$_]) } 0..2;
    my $t;
    for(0..2){
	$t+=$dv[$_];
    }
    $t>1 ? undef : 1;
}

######################################################################
# Mutation strategies
######################################################################
sub end_rotation {
    my $pkg  = shift;
    my $fold = shift;
    my ($node, $end);

    if(int(rand(10000)%2)){
	$node = 1;
	$end = 0;
    }
    else{
	$node = $pkg->{length} - 2;
	$end = $node + 1;
    }

    my $old_end_pos = $fold->{chain}->[$end];
    my $old_node_pos = $fold->{chain}->[$node];

    foreach my $dir (shuffle @dir){
	my $new_pos = one_step($dir, $fold->{chain}->[$node]);

	if( !pos_equal($new_pos, $old_node_pos) && is_pos_free($fold, $new_pos) ){
	    $fold->{chain}->[$end] = \@$new_pos;
	    $fold->{idx}->{join q/,/, @{$new_pos}} =
 $fold->{idx}->{join q/,/, @{$old_end_pos}};
	    delete $fold->{idx}->{join q/,/, @{$old_end_pos}};
	    last;
	}
    }
}

sub can_do_kink_jump {
    my $fold = shift;
    my $nodeidx = shift;
    my $prev = $fold->{chain}->[$nodeidx - 1];
    my $this = $fold->{chain}->[$nodeidx];
    my $next = $fold->{chain}->[$nodeidx + 1];

    my $a = [ map { $next->[$_] - $this->[$_] } 0..2 ];
    my $b = [ map { $prev->[$_] - $this->[$_] } 0..2 ];
    my $t = [ map { $a->[$_] + $b->[$_] } 0..2 ];
    my $t2 =[ map { $this->[$_] + $t->[$_] } 0..2 ];
    (is_pos_free($fold, $t2), $t2);
}

sub kink_jump {
    my $pkg  = shift;
    my $fold = shift;
    for my $n (1..$pkg->{length}-2){
	my ($cando, $pos) = can_do_kink_jump($fold, $n);
	if(int(rand(10000)%2) && $cando){
	    delete $fold->{idx}->{join q/,/, @{$fold->{chain}->[$n]}};
	    $fold->{idx}->{join q/,/, @{$pos}} = $pkg->{typearr}->[$n];
	    $fold->{chain}->[$n] = $pos;
	    last;
	}
    }
}

sub can_do_crankshaft {
    my $fold = shift;
    my $nodeno = shift;
    return undef unless is_dist_one($fold->{chain}->[$nodeno],
				    $fold->{chain}->[$nodeno + 3]); 
    my($tmppos, $tmppos2);
    $tmppos = select_freepos($fold, $fold->{chain}->[$nodeno]);
    if(is_pos_free($fold, $tmppos)){
	$tmppos2 = posadd(
			  possub($tmppos, $fold->{chain}->[$nodeno+1]),
			  $fold->{chain}->[$nodeno+2]
			  );
	if(is_pos_free($fold, $tmppos2)){
	    return (1, $tmppos, $tmppos2);
	}
    }
    undef;
}

sub crankshaft {
    my $pkg  = shift;
    my $fold = shift;

    foreach my $n (0..$pkg->{length}-4){
	my ($cando, $pos1, $pos2) = can_do_crankshaft($fold, $n);
	if(int(rand(10000)%2) && $cando){
	    delete $fold->{idx}->{join q/,/, @{$fold->{chain}->[$n+1]}};
	    delete $fold->{idx}->{join q/,/, @{$fold->{chain}->[$n+2]}};
	    $fold->{idx}->{join q/,/, @{$pos1}} = $pkg->{typearr}->[$n+1];
	    $fold->{idx}->{join q/,/, @{$pos2}} = $pkg->{typearr}->[$n+2];

	    $fold->{chain}->[$n+1] = $pos1;
	    $fold->{chain}->[$n+2] = $pos2;
	    last;
	}

    }
}

sub slithering {
    my $pkg  = shift;
    my $fold = shift;
    my $new_head_pos;
    my $tmppos;

    if(int(rand(10000)%2)){
	$new_head_pos = select_freepos($fold, $fold->{chain}->[$pkg->{length}-1]);
	foreach my $id (reverse 0..$pkg->{length}-1){
	    my $oldpos = $fold->{chain}->[$id];
	    if($id == 0){
		delete $fold->{idx}->{join q/,/, @$oldpos};
	    }

	    $fold->{chain}->[$id] = $id == $pkg->{length}-1 ? $new_head_pos : $tmppos;
	    $fold->{idx}->{join q/,/, @{$fold->{chain}->[$id]}} = $pkg->{typearr}->[$id];
	    $tmppos = $oldpos;
	}
    }
    else {
	$new_head_pos = select_freepos($fold, $fold->{chain}->[0]);
	foreach my $id (0..$pkg->{length}-1){
	    my $oldpos = $fold->{chain}->[$id];
	    if($id == $pkg->{length}-1){
		delete $fold->{idx}->{join q/,/, @$oldpos};
	    }
	    $fold->{chain}->[$id] = $id == 0 ? $new_head_pos : $tmppos;
	    $fold->{idx}->{join q/,/, @{$fold->{chain}->[$id]}} = $pkg->{typearr}->[$id];
	    $tmppos = $oldpos;
	}

    }

}


######################################################################
# Energy Function
######################################################################

sub calc_energy {
    my $pkg = shift;
    my $fold = $pkg->{pool}->[shift];
    my @pos;
    $fold->{energy} = 0;
    foreach (keys %{$fold->{idx}}){
      @pos = split /,/;
      foreach my $i (0..5){
	my $tmp = one_step($dir[$i], \@pos);
	if($fold->{idx}->{join q/,/, @{$tmp}} == 1 &&
	   $fold->{idx}->{join q/,/, @pos} == 1){
	  $fold->{energy} -= 1;
	}
      }
    }
    $fold->{energy} /= 2;
}

######################################################################
# Mutation method
######################################################################

# foreach (qw/end_rotation kink_jump crankshaft slithering/){
#     print <<CODE;
# if(\$r < \$pkg->{mutdist}->{$_}){
#   print "$_\n";
#  \$pkg->$_(\$pkg->{pool}->[\$foldidx]);
#   return;
# }
# CODE
# }


sub mutate_one {
    my $pkg  = shift;
    my $foldidx = shift;
    my $r = rand(1);

    if($r < $pkg->{mutdist}->{end_rotation}){
	$pkg->end_rotation($pkg->{pool}->[$foldidx]);
	return;
    }
    if($r < $pkg->{mutdist}->{kink_jump}){
	$pkg->kink_jump($pkg->{pool}->[$foldidx]);
	return;
    }
    if($r < $pkg->{mutdist}->{crankshaft}){
	$pkg->crankshaft($pkg->{pool}->[$foldidx]);
	return;
    }
    $pkg->slithering($pkg->{pool}->[$foldidx]);
}

sub mutate {
    my $pkg = shift;
    foreach (0..$pkg->{popsize}-1){
	$pkg->mutate_one($_);
	$pkg->calc_energy($_);
    }
}


sub refresh_pool {
  my $pkg = shift;
  my $energy_thresold =
      (sort { $a <=> $b } 
       map { $pkg->{pool}->[$_]->{energy} } 0..$pkg->{popsize}-1)
	  [$pkg->{survival}];

  for (0..$pkg->{popsize}-1){
      $pkg->{pool}->[$_] = undef if $pkg->{pool}->[$_]->{energy} > $energy_thresold;
      $pkg->create_fold($_) unless defined $pkg->{pool}->[$_];
  }
}

sub run {
  my $pkg = shift;
  $pkg->create_pool();
  foreach (1..$pkg->{cycle}){
      $pkg->mutate();
      $pkg->refresh_pool();
  }

}

sub result {
    my $pkg = shift;
    map { { chain => $_->{chain}, energy => $_->{energy} } }
    grep { defined $_->{energy} }
    @{$pkg->{pool}};
}

1;
__END__
# Below is stub documentation for your module. You better edit it!

=head1 NAME

Bio::Tools::LatticeProtein - Lattice model for protein prediction

=head1 SYNOPSIS

  use Data::Dumper;
  use Bio::Tools::LatticeProtein;
  $lp = Bio::Tools::LatticeProtein->new(
					sequence       => 'PERLISNIFTY',
				        cycle          => 100,
                                        popsize        => 20,
                                        survival       => 10,
					mutdist        => {
					    end_rotation => 0.1,
					    kink_jump    => 0.4,
					    crankshaft   => 0.4,
					    slithering   => 0.1,
					},
					);
  print $lp->run;
  print Dumper $lp->result;


=head1 DESCRIPTION

This module is an experimental implementation of the lattice model for protein tertiary structure prediction. Positions of amino acids are mapped to discrete coordinates. With the genetic algorithm backend, structures are mutated and those with lower energies will be picked up for next round of mutation.

There are four supported mutation strategies. I<End Rotation>, I<Kink Jump>, I<Crankshaft>, and I<Slithering Snake>.


Then, follow the 2D-graphic illustrations of these operations.

=head2 End Rotation


                                  o
                                  |
                =========>        |
                                  |
 o-----o-----o                    o-----o

=head2 Kink Jump		
						
						
  o				   o-----o	
  |    	       	       	       	    	 |	
  |             =========>	    	 |	
  |				    	 |	
  o-----o-----o                          o-----o


=head2 Crankshaft

                                   o-----o
                                   |     |
                                   |     |
                =========>         |     |
  o     o-----o                    o     o-----o
   \     \
    \     \
     \     \
      o-----o

=head2 Slithering snake

  o-----o-----o                                  o-----o
              |                                        |
              |                                        |
              |                   ==========>          |
              o     O <- The head                      o     o-----O
              |     |                                  |     |
              |     |                                  |     |
              |     |                                  |     |
              o-----o                                  o-----o


=head1 INTERFACE

=head2 new

The constructor takes the following parameters

=head3 sequence

The protein sequence. Case-insensitive.

=head3 cycle

The number of generations will be performed. Defaults to 15.

=head3 popsize

The population size of a generation. Defaults to 20.

=head3 survival

The number of protein folding configurations will be carried to next cycle. Defaults to 10, and it should not be greater than popsize.

=head3 mutdist

The distribution of mutation probabilities. The Default setting is 

    mutdist        => {
	end_rotation => 0.25,
	kink_jump    => 0.25,
	crankshaft   => 0.25,
	slithering   => 0.25,
    },

And they must sum to 1.


=head2 run

Starts the whole process.

=head2 result

Returns the predicted result with coordinates and energies.

=head1 ON ENERGY

The method to calculate energy is: First, setting it to zero. For each hydrophobic amino acid, if there's another hydrophobic amino acid found around it, then the energy decreases by 1. 

=head1 REFERENCE

I<MOLECULAR MODELLING, principles and applications> by B<Andrew R. Leach>.


=head1 COPYRIGHT

xern E<lt>xern@cpan.orgE<gt>

This module is free software; you can redistribute it or modify it under the same terms as Perl itself.

=cut
