use ExtUtils::MakeMaker;


WriteMakefile(
    'NAME'		=> 'Bio::Tools::LatticeProtein',
    'VERSION_FROM'	=> 'LatticeProtein.pm', 
    'PREREQ_PM'		=> {
	'List::Util' => '1',
    }, 
    ($] >= 5.005 ?    
      (ABSTRACT_FROM => 'LatticeProtein.pm', 
       AUTHOR     => 'xern <xern@cpan.org>') : ()),
);
