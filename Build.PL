use 5.014;
use strict;
use warnings FATAL => 'all';
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'RNAhelper',
    license             => 'gpl_3',
    dist_author         => q{Felix Kuehnl <felix@bioinf.uni-leipzig.de>},
    dist_version_from   => 'lib/RNAhelper.pm',
    release_status      => 'stable',
    configure_requires => {
        'Module::Build' => '0',
    },
    build_requires => {
        'Test::More' => '0',
    },
    requires => {
        #'ABC'              => '1.6',
        #'Foo::Bar::Module' => '5.0401',
    },
    add_to_cleanup     => [ 'RNAhelper-*' ],
);

$builder->create_build_script();
