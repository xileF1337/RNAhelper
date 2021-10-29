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
        'Carp'                   => '1.38',
        'warnings'               => '1.37',
        'autodie'                => '2.29',
        'IPC::System::Simple'    => '1.25',
        'experimental'           => '0.016',
        'Exporter'               => '5.72',
        'Readonly'               => '2.05',
        'Scalar::Util'           => '1.5',
        'List::Util'             => '1.5',
        'List::MoreUtils'        => '0.428',
        'Math::Round'            => '0.07',
        'Clone'                  => '0.41',
        'File::Spec'             => '3.75',
        'File::Basename'         => '2.85',
        'RNA'                    => '2.4.13', # NOT on cpan, install manually
        'Devel::Assert'          => '1.06',
        'Math::Random::MT::Auto' => '6.22',
        'Inline'                 => '0.83',
        'Inline::C'              => '0.78',
        'Math::BigFloat'         => '1.999806',
    },
    add_to_cleanup     => [ 'RNAhelper-*' ],
);

$builder->create_build_script();
