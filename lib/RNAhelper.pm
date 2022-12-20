# RNAhelper.pm
# Version 1.4 2022-12-20
#
# Provide utilities like pairtable generation, neighbor generation,
# gradient walks etc. for RNA secondary structures.
#
# Perform a gradient walk among all neighbor RNA structures.
# Move set: single base pair addition / removal, shift moves, lonely pair
#   removal
# Requires ViennaRNA Perl bindings
#
#
# Copyright 2015--2022 Felix Kuehnl, felix[at]bioinf.uni-leipzig.de
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################


package RNAhelper;

use v5.12;
use warnings;
use autodie ':all';

our $VERSION = '1.4';

#use diagnostics;               # REALLY verbose error msges
# Disables warning independent of perl version, requires module 'experimental'
use experimental 'smartmatch';

use Exporter;

# use Data::Dumper;
use Carp qw( croak confess );
use Readonly;
use Scalar::Util qw( blessed reftype );
use List::Util qw( max min reduce pairs none );
use List::MoreUtils qw( any );
use Math::Round;
# Shit, Clone::Fast seems to be broken with Perl >5.22
# use Clone::Fast qw( clone );      # deep copy perl data structures, faster
use Clone 'clone';
# Priority queue, call $ref=new Heap:Priority, ->pop/add, WARNING: modified version, added peek()!
# Use Mersenne twister for pseudo-random number generation.
use Math::Random::MT::Auto qw(rand srand irand get_seed);

use File::Spec;                 # catfile()
use File::Basename;             # dirname

use RNA;                        # ViennaRNA Perl bindings
use Devel::Assert 'off';        # use 'on/'off' to switch assertions on/off

# Install modules Inline and Inline::C!


########################
##  Global variables  ##
########################

# Add all valid BPs. Add both xy and yx for symmetry.
# use Readonly;
# Readonly::Array my @VALID_BP => ('AU','UA','GC','CG','GU','UG');
my @VALID_BP = qw(AU UA GC CG GU UG);

# Universal gas constant in kcal/mol as in ViennaRNA
# my $GASCONST = 1.98717e-3;
Readonly::Scalar my $GASCONST => RNA::GASCONST * 1e-3;


#####################
##  Package setup  ##
#####################

BEGIN {
    ### Module options
    our @ISA=('Exporter');
    # These subs WILL be exported if module is loaded
    our @EXPORT = qw();
    # These subs CAN be exported if explicitly specified
    our @EXPORT_OK = qw(
        printStrNrg
        printStrNrgLst
        pt2str pt2str_c struct_string
        str2pt
        path2str

        rand_seq
        rand_seed

        grad_walk
        floodCycle

        bestDirectPath
        bestDirectSaddle
        directPaths

        genNeighbors
        genNeighbors_pairs
        isLonely            is_lonely
        isLonelyStruct      is_lonely_struct
        isCanonical         is_canonical
        isInsLonely         is_ins_lonely
        isOutLonely         is_out_lonely

        is_valid_seq
        normalize_seq
        seq_len

        bpValid         bp_valid
        is_paired
        struct_cmp
        struct_lt
        struct_eq
        bp_dist
        bp_iter_factory
        nrg_of_move

        rna_mfe
        rna_mfe_struct
        rna_eval_structure
        make_fold_compound
        make_model_details

        boltz2nrg
        nrg2boltz
        ensemble_energy
        subopt_energy
        partition_energy
        partition_energy_structiter
        partition_energy_nrgiter
        get_gas_const
        celsius2kelvin
        get_basin_partition_funcs

        base_pair_prob

        get_opt
        parUniq
        parSort

        test
    );
}


#assert test();
test() unless caller;       # Test only if module is executed directly


#############
##  Tests  ##
#############

# Convert dot bracket strings from test input arrays to pairtable references
sub convTstInput( @ )
{
    my @ret = @_;
    for( my $i=0; $i<@ret; $i+=2)
    {
        $ret[$i] = str2pt($ret[$i]);
    }
    return @ret;
}


# Round energies in struct/nrg lists to two digits after decimal point
sub roundNrgs( @ )
{
    my @ret = @_;
    for( my $i=1; $i<@ret; $i+=2)
    {
        $ret[$i] = nearest( 0.01, $ret[$i]);
    }
    return @ret if wantarray;
    return \@ret;
}


sub cmpMoves
{
    my ($a, $b) = @_;
    if( @$a < @$b) { return -1 }
    elsif (@$a > @$b) { return 1 }
    for my $i (0..$#$a) {
        if(    $a->[$i][0] < $b->[$i][0]) { return -1 }
        elsif( $a->[$i][0] > $b->[$i][0]) { return  1 }
        elsif( $a->[$i][1] < $b->[$i][1]) { return -1 }
        elsif( $a->[$i][1] > $b->[$i][1]) { return  1 }
    }
    return 0;
}

# To run this test suite, simply execute 'perl Module_name.pm'.
sub test
{
    # Load testing framework, skip_all => "reason' or to skip
    eval { require Test::More }
        or croak "Module Test::More is required to run the test suite";
    Test::More->import();

    print "Hello from the test suite, wise are you to run me :-)\n";


    # isInsLonely, isOutLonely, isLonely, outGrowLonely, insGrowLonely, isLonelyStruct
    my $seq = 'GGGGUUUCCCC';
    my @str = qw{
        .(((....))).
        ..((....))..
        ..(......)..
        .((......).)
        .((.(...)).)
        .(((..()))).
        (((......)))
        ((((.)(.))))
    };
    #   012345678901

    my @insIsLone      = (0, 0, 1, 1, 1, 0, 1, 1);
    my @outIsLone      = (0, 1, 1, 1, 1, 0, 0, 0);
    my @isLone         = (0, 0, 1, 1, 1, 0, 0, 0);
    my @insGrowLone    = (1, 1, 0, 0, 0, 1, 0, 0);
    my @outGrowLone    = (1, 0, 0, 0, 0, 1, 0, 0);
    my @isLonelyStruct = (0, 0, 1, 1, 1, 0, 0, 1);

    ok( @str==@insIsLone && @str==@outIsLone && @str==@insGrowLone
               && @str==@outGrowLone && @str==@isLone,
           'loneliness tests: array lengths');
    for my $i( 0..$#str)
    {
        my @pt = str2pt( $str[$i]);
        cmp_ok( isInsLonely(   \@pt, 2, 9), '==', $insIsLone[$i],   "isInsLonely(2,9)   - $str[$i]");
        cmp_ok( isOutLonely(   \@pt, 2, 9), '==', $outIsLone[$i],   "isOutLonely(2,9)   - $str[$i]");
        cmp_ok( isLonely(      \@pt, 2, 9), '==', $isLone[$i],      "isLonely(2,9)      - $str[$i]");
        cmp_ok( insGrowLonely( \@pt, 2, 9), '==', $insGrowLone[$i], "insGrowLonely(2,9) - $str[$i]");
        cmp_ok( outGrowLonely( \@pt, 2, 9), '==', $outGrowLone[$i], "outGrowLonely(2,9) - $str[$i]");
        cmp_ok( isLonelyStruct( $seq, $str[$i]), '==', $isLonelyStruct[$i],
                                                                    "isLonelyStruct     - $str[$i]");
    }

    ### genNeighbors
    my $str = '((((...))))';
    my @pt  = str2pt( $str);
    # Correct restoration of symbol table when using genNeighbors with
    # do_energies option set to false.
    my @nb  = roundNrgs genNeighbors( $seq, \@pt, 'Inf', 0, 1, undef, undef, 0);
    # Actual neighbor computation
    @nb  = roundNrgs genNeighbors( $seq, \@pt, 9999, 0, 0);
    my @res =  convTstInput (
        '(((.....)))', -2.40,
        '.(((...))).', -2.20,
        '(.((...)).)',  2.60,
        '((.(...).))',  2.60,
       );
    is_deeply( \@nb, \@res, 'genNeigbors: rem noLP=0 shift=0');

    @nb  = roundNrgs genNeighbors( $seq, \@pt, 9999, 1, 0);
    @res =  convTstInput (
        '(((.....)))', -2.40,
        '.(((...))).', -2.20,
       );
    is_deeply( \@nb, \@res, 'genNeigbors: rem noLP=1 shift=0');

    #seq = 'GGGGUUUCCCC';
    $str = '...........';
    @pt  = str2pt( $str);
    @nb  = roundNrgs genNeighbors( $seq, \@pt, 9999, 1, 0);    # minh=inf, noLP, no slide moves
    @res =  convTstInput (
           '((......)).', -0.20, '.((....))..', -0.20, '.((.....)).', -0.10,
           '((.....))..',  0.10, '..((....)).',  0.30, '.((......))',  0.60,
           '..((...))..',  1.10, '((.......))',  1.20, '..((.....))',  1.40,
           '((....))...',  2.80, '.((...))...',  3.40, '((...))....',  5.80,
       );
    is_deeply( \@nb, \@res, 'genNeigbors: add noLP=1 shift=0');

    $str = '..((...))..';
    @pt  = str2pt( $str);
    @nb  = roundNrgs genNeighbors( $seq, \@pt, 9999, 1, 0);    # minh=inf, noLP, no slide moves, bug 1
    @res =  convTstInput (
        '.(((...))).', -2.20,
        '...........',  0.00,
    );
    is_deeply( \@nb, \@res, 'genNeigbors: add noLP=1 shift=0, bug 1');

    $seq = 'GGGGGGGGGGGGAAAUUUUUUUGUUUU';
    $str = '.(((....(((.....)))....))).';
    #       01234567890123456789012
    #       0         1         2
    @pt  = str2pt( $str);
    @nb  = roundNrgs genNeighbors( $seq, \@pt, 9999, 0, 1);    # minh=inf, noLP=0, shift moves
    @res =  convTstInput(
        '.(((...((((.....))))...))).',  4.90, '.(((.(..(((.....)))..).))).',  5.40,
        '((((....(((.....)))....))))',  5.90, '.(((....((((...))))....))).',  5.90,
        '..((....(((.....)))....))..',  6.20, '.(((.....((.....)).....))).',  6.40,
        '.(((....((.......))....))).',  6.50, '.(((..(.(((.....))).)..))).',  6.70,
        '.((.....(((.....))).....)).',  7.40, '.(.(....(((.....)))....).).',  8.60,
        '.(((....(.(.....).)....))).',  8.60, '.(((....(((....).))....))).',  8.90,
        '.(((.(..(((.....))).)..))).',  9.30, '.(((..(.(((.....)))..).))).',  9.30,
        '.(((....((.(....)))....))).',  9.40, '.(((..(.(((.....))))...))).',  9.60,
        '.(((...((((.....))).)..))).',  9.60, '.(((....(((.....)))....)).)',  9.80,
        '.(((...(.((.....)))....))).',  9.90, '.(((....(((.....)).)...))).',  9.90,
        '.((.(...(((.....)))....))).',  9.90, '(.((....(((.....)))....))).', 10.10,
        '.((((...(((.....)))..).))).', 10.30, '.((((...(((.....))).)..))).', 10.50,
        '.(((.(..(((.....))))...))).', 10.50, '.(((...((((.....)))..).))).', 10.50,
        '.(((....((..(...)))....))).', 10.70, '.(((..(..((.....)))....))).', 10.90,
        '.(((....(((.....))..)..))).', 10.90, '.((..(..(((.....)))....))).', 10.90,
        '.((((...(((.....))))...))).', 11.40, '.(((....(((.....)))..)..)).', 11.90,
        '.((((....((.....)))....))).', 12.70, '.((....((((.....)))....))).', 12.70,
        '.(((....(((.....))))....)).', 12.70, '.(((.(...((.....)))....))).', 12.90,
        '.(((....(((.....))...).))).', 12.90, '.((...(.(((.....)))....))).', 12.90,
        '.(((....(((.....))).)...)).', 12.90, '.(((.....((.....))(...)))).', 15.40,
    );
    is_deeply( \@nb, \@res, 'genNeigbors: shift noLP=0 shift=1');

    $seq = 'GGGGGGGGUUUUUUU';
    $str = '((..(((...)))))';
    #       01234567890123456789012
    #       0         1         2
    @pt  = str2pt( $str);
    @nb  = roundNrgs genNeighbors( $seq, \@pt, 9999, 1, 1);    # minh=inf, noLP=1, shift moves
    @res =  convTstInput (
           '....(((...)))..', 4.70,
           '((..((.....))))', 8.50,
           '(((..((...)))))', 8.70,
           '((...((...)).))', 9.10,
    );
    is_deeply( \@nb, \@res, 'genNeigbors: shift noLP=1 shift=1');

    $seq = 'AAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUUACAUCUGAAGUGCUGCCUUUUUUUU';
    $str = '(((((((......)))).)))((((((...)))).....)).........((((......)))).....';
    @pt  = str2pt( $str);
    @nb  = roundNrgs genNeighbors( $seq, \@pt, 0.5, 1, 1);    # minh=inf, noLP=1, shift moves
    @res =  convTstInput (
        '(((((((......)))).)))..((((...))))................((((......)))).....', -5.3,
        '.((((((......)))).))(((((((...)))).....)))........((((......)))).....', -3.9
    );
    is_deeply( \@nb, \@res, 'genNeigbors: shift noLP=1 shift=1 cross shifts');

    $seq = 'AAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAU';
    $str = '(((((((......)))).)))((((((...)))).....))..';
    @pt  = str2pt( $str);
    @nb  = roundNrgs genNeighbors( $seq, \@pt, 9999, 1, 1);    # minh=inf, noLP=1, shift moves
    @res =  convTstInput (
        '(((((((......)))).)))..((((...)))).........',  -5.40,
        '.((((((......)))).))(((((((...)))).....))).',  -4.00,
        '.((((((......)))).)).((((((...)))).....))..',  -3.10,
        '((((((........))).)))((((((...)))).....))..',  -2.60,
        '(((.(((......)))..)))((((((...)))).....))..',  -2.40,
        '(((((((......)))).)))(((((.....))).....))..',  -1.40,
        '(((((((......)))).)))((.(((...)))......))..',  -1.10,
        '((.((((......))))..))((((((...)))).....))..',   0.50,
        '(((((((......))))).))((((((...)))).....))..',   1.30,
    );
    is_deeply( \@nb, \@res, 'genNeigbors: shift noLP=1 shift=1 cross shifts2');

    $seq = 'AAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUUUACAUACUCG';
    $str = '(((((...((((....((......))..))))...)))))..............';
    @pt  = str2pt( $str);
    @nb  = roundNrgs genNeighbors( $seq, \@pt, 0, 0, 1);    # minh=0, noLP=0, shift moves
    @res =  convTstInput (
        # This can be empty, as the 'evil structure'
        #   (((((...(((((...)(......)...))))...)))))..............   4.30
        # is generated anyway, but not reported.
    );
    is_deeply( \@nb, \@res, 'genNeigbors: shift noLP=0 shift=1 bug: ViennaRNA dies if j>i (cross shift moves)');


    # say "Printing results\n", $str;
    # while( my @nxt = splice( @nb, 0, 2))
    # {
    #     say printStrNrg( $nxt[0], $nxt[1]);
    # }

    # for( my $i=0; $i<@nb; $i++)
    # {
    #     local $, = ' ';
    #     say "-$i-";
    #     say "nb:  ", $i%2 ? $nb[$i] : @{$nb[$i]};
    #     say "res: ", $i%2 ? $res[$i] : @{$res[$i]};
    # }


    ################## ## Direct paths ## ##################

    $seq      = "GGGUUUUUCCC";
    my $stra  = "(((...))..)";
    my $strb  = "((.......))";
    #            01234567890
    my @paths = sort { cmpMoves($a->{moves}, $b->{moves}) }
                       directPaths( $seq, [str2pt( $stra)], [str2pt( $strb)], 50,
                                       { shift=>0, noLP=>0 });

    @res = (
        [ [1,7], [1,9], [2,6], ],
        [ [1,7], [2,6], [1,9], ],
        [ [2,6], [1,7], [1,9], ],
    );
    is_deeply( [ map { $_->{moves} } @paths ], \@res, 'directPaths: noLP=0 shift=0, no crossing');

    @paths = sort { cmpMoves($a->{moves}, $b->{moves}) }
                directPaths( $seq, [str2pt( $stra)], [str2pt( $strb)], 50,
                                   { shift=>1, noLP=>0 });
    @res = (
        [ [1,7], [1,9], [2,6], ],
        [ [1,7], [1,9], [2,6], ],
        [ [1,7], [2,6], [1,9], ],
        [ [2,6], [1,7], [1,9], ],
        [ [2,6], [1,7], [1,9], ],
    );
    is_deeply( [ map { $_->{moves} } @paths ], \@res, 'directPaths: noLP=0 shift=1, no crossing');

    $seq  = "GGGUUUUCCC";
    $stra = ".((.....))";
    $strb = "((....))..";
    #        01234567890
    @paths = sort { cmpMoves($a->{moves}, $b->{moves}) }
             directPaths( $seq, [str2pt( $stra)], [str2pt( $strb)], 50,
                               { shift=>0, noLP=>0 });
    @res = (
        [ [1,9], [2,8], [0,7], [1,6], ],
        [ [1,9], [2,8], [1,6], [0,7], ],
        [ [2,8], [1,9], [0,7], [1,6], ],
        [ [2,8], [1,9], [1,6], [0,7], ],
    );
    is_deeply( [ map { $_->{moves} } @paths ], \@res, 'directPaths: noLP=0 shift=0, crossing');

    @paths = sort { cmpMoves($a->{moves}, $b->{moves}) }
             directPaths( $seq, [str2pt( $stra)], [str2pt( $strb)], 50,
                               { shift=>1, noLP=>0 });
    @res = (
        [ [1,9], [2,8], [0,7], [1,6], ],
        [ [1,9], [2,8], [1,6], [0,7], ],
        [ [2,8], [1,9], [0,7], [1,6], ],
        [ [2,8], [1,9], [1,6], [0,7], ], # yea they are equal, two middle bp are a double move
        [ [2,8], [1,9], [1,6], [0,7], ],
    );
    is_deeply( [ map { $_->{moves} } @paths ], \@res, 'directPaths: noLP=0 shift=1');

    @paths = sort { cmpMoves($a->{moves}, $b->{moves}) }
             directPaths( $seq, [str2pt( $stra)], [str2pt( $strb)], 50,
                               { shift=>0, noLP=>1 });
    @res = (
        [ [1,9], [2,8], [0,7], [1,6], ],
    );
    is_deeply( [ map { $_->{moves} } @paths ], \@res, 'directPaths: noLP=1 shift=0, crossing');

    $seq =  'GGGGGGGGGGGGGGGCCCCCCCCCCCCCC';
    $stra = '.(((....(((.....))).....)))..';
    $strb = '.((......((((((....)))).)).))';
    #        01234567890
    @paths = sort { cmpMoves($a->{moves}, $b->{moves}) }
             directPaths( $seq, [str2pt( $stra)], [str2pt( $strb)], 100,
                               { shift=>0, noLP=>1, nreturn=>10 });
    @res = (
        [
            [10,16], [8,18], [9,17], [11,22], [12,21], [13,20], [ 1,26],
            [ 2,25], [3,24], [9,25], [10,24], [ 1,28], [ 2,27], [14,19],
        ],
        [
            [10,16], [8,18], [9,17], [11,22], [12,21], [13,20], [1,26],
            [2,25], [3,24], [9,25], [10,24], [14,19], [1,28], [2,27],
        ],
        [
            [10,16], [8,18], [9,17], [11,22], [12,21], [13,20], [1,26],
            [2,25], [3,24], [14,19], [1,28], [2,27], [9,25], [10,24],
        ],
        [
            [10,16], [8,18], [9,17], [11,22], [12,21], [13,20], [1,26],
            [2,25], [3,24], [14,19], [9,25], [10,24], [1,28], [2,27],
        ],
        [
            [10,16], [8,18], [9,17], [11,22], [12,21], [13,20], [1,26],
            [14,19], [2,25], [3,24], [1,28], [2,27], [9,25], [10,24],
        ],
        [
            [10,16], [8,18], [9,17], [11,22], [12,21], [13,20], [1,26],
            [14,19], [2,25], [3,24], [9,25], [10,24], [1,28], [2,27],
        ],
        [
            [10,16], [8,18], [9,17], [11,22], [12,21], [13,20], [14,19],
            [1,26], [2,25], [3,24], [1,28], [2,27], [9,25], [10,24],
        ],
        [
            [10,16], [8,18], [9,17], [11,22], [12,21], [13,20], [14,19],
            [1,26], [2,25], [3,24], [9,25], [10,24], [1,28], [2,27],
        ],
        [
            [10,16], [8,18], [9,17], [11,22], [12,21], [13,20], [14,19],
            [3,24], [1,26], [2,25], [1,28], [2,27], [9,25], [10,24],
        ],
        [
            [10,16], [8,18], [9,17], [11,22], [12,21], [13,20], [14,19],
            [3,24], [1,26], [2,25], [9,25], [10,24], [1,28], [2,27],
        ],
    );
    is_deeply( [ map { $_->{moves} } @paths ], \@res, 'directPaths: noLP=1 shift=0, no. 2');

    $seq =  'GGGGGGGGGCCCCCCCCCCC';
    $stra = '....((((.......)))).';
    $strb = '....(((((....))).)).';
    #        01234567890
    @paths = sort { cmpMoves($a->{moves}, $b->{moves}) }
             directPaths( $seq, [str2pt( $stra)], [str2pt( $strb)], 100,
                               { shift=>1, noLP=>0, nreturn=>20});
    @res = (
        [ [6,16], [7,15], [6,15], [7,14], [8,13], ], [ [6,16], [7,15], [6,15], [7,14], [8,13], ],
        [ [6,16], [7,15], [6,15], [8,13], [7,14], ], [ [6,16], [7,15], [6,15], [8,13], [7,14], ],
        [ [7,15], [6,16], [6,15], [7,14], [8,13], ], [ [7,15], [6,16], [6,15], [7,14], [8,13], ],
        [ [7,15], [7,14], [6,16], [6,15], [8,13], ], [ [7,15], [7,14], [6,16], [6,15], [8,13], ],
        [ [7,15], [7,14], [8,13], [6,16], [6,15], ], [ [7,15], [7,14], [8,13], [6,16], [6,15], ],
        [ [7,15], [7,14], [8,13], [6,16], [6,15], ], [ [7,15], [7,14], [8,13], [6,16], [6,15], ],
        [ [7,15], [8,13], [6,16], [6,15], [7,14], ], [ [7,15], [8,13], [7,14], [6,16], [6,15], ],
        [ [7,15], [8,13], [7,14], [6,16], [6,15], ], [ [8,13], [7,15], [6,16], [6,15], [7,14], ],
        [ [8,13], [7,15], [7,14], [6,16], [6,15], ], [ [8,13], [7,15], [7,14], [6,16], [6,15], ],
        [ [8,13], [7,15], [7,14], [6,16], [6,15], ], [ [8,13], [7,15], [7,14], [6,16], [6,15], ],
    );
    is_deeply( [ map { $_->{moves} } @paths ], \@res, 'directPaths: noLP=0 shift=1, no. 3');

    $seq =  'AAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUUACAUCUGAAGUGCUGCCUUUUUUUU';
    $stra = '.(((((.....((..........((((...))))))...))))).........................';
    $strb = '(((((((......)))).)))..((((...))))......((.(((.......))).))..........';
    #        0123456789012345678901234567890123456789012345678901234567890123456789
    #        0         1         2         3         4         5         6         7
    @paths = path2str( (directPaths( $seq, [str2pt( $stra)], [str2pt( $strb)], 30,
                               { shift=>0, noLP=>0, nreturn=>1}))[0]);
    @res = (
        '.(((((.....((..........((((...))))))...))))).........................  -4.80',
        '.(((((......(..........((((...)))))....))))).........................  -1.20',
        '..((((......(..........((((...)))))....))))..........................  -0.30',
        '..((((.................((((...)))).....))))..........................  -3.30',
        '..((((.................((((...)))).....))))..(.......)...............   0.10',
        '..((((.................((((...)))).....)))).((.......))..............  -1.30',
        '...(((.................((((...)))).....)))..((.......))..............   0.20',
        '...(((.................((((...)))).....))).(((.......))).............  -0.70',
        '....((.................((((...)))).....))..(((.......))).............   0.00',
        '....(..................((((...))))......)..(((.......))).............   2.40',
        '.......................((((...)))).........(((.......))).............  -2.60',
        '....(..........).......((((...)))).........(((.......))).............   0.90',
        '....(.(......).).......((((...)))).........(((.......))).............   2.20',
        '....(((......))).......((((...)))).........(((.......))).............  -2.50',
        '..(.(((......)))..)....((((...)))).........(((.......))).............  -0.50',
        '.((.(((......)))..))...((((...)))).........(((.......))).............  -1.80',
        '.((.(((......)))..))...((((...))))......(..(((.......)))..)..........   0.60',
        '(((.(((......)))..)))..((((...))))......(..(((.......)))..)..........  -0.30',
        '(((.(((......)))..)))..((((...))))......((.(((.......))).))..........  -3.20',
        '(((((((......)))).)))..((((...))))......((.(((.......))).))..........  -4.80',
    );
    is_deeply( [ @paths ], \@res, 'directPaths: bug: not all required bp added to add list');

    #$RNA::temperature = 37;
    # say "$stra\n$strb\n01234567890123456789";
    # for my $path (@paths)  # Print paths
    # {
    #     #print "Path ", $i++, ": ";
    #     print "[";
    #     print " [$_->[0],$_->[1]]," foreach @{$path->{moves}};
    #     print " ],\n";
    #     printf "saddle: %6.2f\n", $path->{saddle};
    #     say map { sprintf "%6.2f ", $_ } @{$path->{nrgs}};
    #     say for path2str( $path);
    # }

    $seq =  'GGGGGGGGGCCCCCCCCCCC';
    $stra = '....((((.......)))).';
    $strb = '....(((((....))).)).';
    #        01234567890
    @paths = sort { cmpMoves($a->{moves}, $b->{moves}) }
             directPaths( $seq, [str2pt( $stra)], [str2pt( $strb)], 100,
                               { shift=>1, noLP=>1, nreturn=>20});
    @res = (
        [ [7,15], [6,16], [6,15], [7,14], [8,13], ],
        [ [7,15], [6,16], [7,14], [8,13], [6,15], ],
        [ [7,15], [7,14], [8,13], [6,16], [6,15], ],
        [ [7,15], [7,14], [8,13], [6,16], [6,15], ],
    );
    is_deeply( [ map { $_->{moves} } @paths ], \@res, 'directPaths: noLP=1 shift=1, no. 3');

    @res = (
        "....((((.......)))).  -6.40",
        "....(((.........))).  -2.70",
        "....(((((....)).))).  -5.50",
        "....(((((....))).)).  -5.50",
    );
    is_deeply( [bestDirectPath( $seq, $stra, $strb, 100, 1, 1)], \@res, "bestDirectPath");
    is( nearest( 0.01, bestDirectSaddle( $seq, $stra, $strb, 100, 1, 1)), -2.7,
        "bestDirectSaddle");

    ############ Partition functions and rates ##########

    my $test_data_dir = 'RNAhelper_test';

    die "test data directory '$test_data_dir' is no valid directory"
         unless -d $test_data_dir;

    my $mfe = -27.4;
    my $test_file = 'get_basin_partition_funcs.RS1.subopt_e25_T37.bar_minh1'
                     . '.nd.sbmap.basin50-99';
    open my $sbmap_handle, '<', "$test_data_dir/$test_file"
            or die "could not open file '$test_file' in directory ",
                   "'$test_data_dir'";
    my @Z = get_basin_partition_funcs( $sbmap_handle, { mfe => $mfe });
    close $sbmap_handle;
    # @Z[ 0 .. 49 ] = (0) x 50;  # Set first basins to 0
    @Z = map { defined $_ ? sprintf( "%10.5g", $_) : ' 'x9 . 0 } @Z;
    @res = map { sprintf "%10.5g", $_ } (
                 0,             0,             0,             0,
                 0,             0,             0,             0,
                 0,             0,             0,             0,
                 0,             0,             0,             0,
                 0,             0,             0,             0,
                 0,             0,             0,             0,
                 0,             0,             0,             0,
                 0,             0,             0,             0,
                 0,             0,             0,             0,
                 0,             0,             0,             0,
                 0,             0,             0,             0,
                 0,             0,             0,             0,
                 0,             0,    # Basin 0 - 49 are zero
          5.31e-15,    3.5402e-15,    2.2388e-15,    4.5895e-15,
        4.5899e-15,    3.6013e-15,    2.4056e-15,    7.1105e-15,
        3.3381e-15,    2.1981e-15,    2.5192e-15,    1.4763e-15,
        1.5847e-15,    2.4094e-15,    1.6021e-15,    1.1638e-15,
        1.9895e-15,    1.2065e-15,    1.4618e-15,    1.1364e-15,
        8.7699e-16,    5.0752e-16,    1.2794e-15,    9.5666e-16,
        7.4024e-16,    3.5154e-16,    5.6824e-16,    3.5817e-16,
        3.2697e-16,    4.1152e-16,    8.0666e-16,    6.1405e-16,
        1.2415e-15,    2.2349e-16,    7.5083e-16,    3.3909e-16,
         2.073e-16,    5.3454e-16,    4.0495e-16,    6.5374e-16,
        5.6845e-16,    2.8091e-16,    4.9214e-16,    3.2319e-16,
        5.0126e-16,    2.3883e-16,    4.3743e-16,    1.4627e-16,
        3.7454e-16,    3.1572e-16,    # Basin 50 - 99, rest omitted
    );
    is_deeply( \@Z, \@res, 'get_basin_partition_funcs');

    ####### Structure comparisons: struct_lt, struct_cmp, struct_eq #######

    # For order '(' < ')' < '.' (i.e. lexicographical / ascii order)
    @str = qw{
        ((...)).....
        .(((....))).
        ..((....))..
        ..(......)..
    };

    # # # For order '.' < '(' < ')'
    # @str = qw{
    #     ..(......)..
    #     ..((....))..
    #     .(((....))).
    #     ((...)).....
    # };


    # result of struct[i] cmp struct[j]
    my $res_struct_cmp = [
        [ 0, -1, -1, -1],
        [ 1,  0, -1, -1],
        [ 1,  1,  0, -1],
        [ 1,  1,  1,  0],
    ];

    # [i][j] == is ith struct lt jth struct?
    my $is_struct_lt = [
        [0, 1, 1, 1],
        [0, 0, 1, 1],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
    ];

    # [i][j] == is ith struct eq jth struct?
    my $is_struct_eq = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ];

    ok( @str==@$res_struct_cmp && @str==@$is_struct_lt
            && @str==@$is_struct_eq,
        'struct comparisons: array lengths');

    my @pts;
    push @pts, [str2pt($_)] foreach @str;
    for my $i ( 0..$#str) {
        for my $j ( 0..$#str) {
            my $other_pt = [ @{ $pts[$j] } ];  # don't check on same array
            cmp_ok( struct_cmp($pts[$i], $other_pt),
                    '==',
                    $res_struct_cmp->[$i][$j],
                    "struct_cmp $i,$j: $str[$i] <=> $str[$j]",
                  );
            cmp_ok( struct_lt($pts[$i], $other_pt),
                    '==',
                    $is_struct_lt->[$i][$j],
                    "struct_lt $i,$j:  $str[$i] <   $str[$j]",
                  );
            cmp_ok( struct_eq($pts[$i], $other_pt),
                    '==',
                    $is_struct_eq->[$i][$j],
                    "struct_eq $i,$j:  $str[$i] ==  $str[$j]",
                  );
        }
    }


    done_testing();
    return 1;
}


#############################
##  Submodule definitions  ##
#############################

# Use fast C implementation for most-used functions. Let the directory
# containing the compiled code reside in the subdirectory '.Inline' in the
# directory of this module.
{
    my $c_code;
    my $bin_store_dir;

    BEGIN {
        $c_code = <<'END_OF_C_CODE';

        /*
         * C version of the heavily used pt2str function which converts a
         * structure's pairtable to a dot-bracket string representation.
         * ... - list of integers representing the pair table
         * returns dot-bracket string (yes, the return type has to be void)
         */
        void pt2str_c( int fst,  ...)
        {
            Inline_Stack_Vars;  /* Prepare stack access */
            int i, j;
            int size = Inline_Stack_Items; /* no. of arguments */
            int *pt = malloc( sizeof(int)*(size));    /* pair table*/
            char *str = malloc( sizeof(char)*(size+1)); /* dot-bracket string */

            /* store pairtable */
            for (i = 0; i < size; i++)
                /* convert SV* stack item to scalar + integer value */
                pt[i] =  SvIV(Inline_Stack_Item(i));

            for( i=0; i<size; i++)
            {
                j = pt[i];
                if( j<0) { str[i] = '.'; }
                else if( i<j) { str[i] = '('; }
                else { str[i] = ')'; }
            }
            str[size] = 0;    /* terminate string with zero */

            /* Return dot-bracket string */
            Inline_Stack_Reset;
            /* SV - scalar value, PV - (char) pointer value == string */
            Inline_Stack_Push(sv_2mortal(newSVpv( str, size)));
            Inline_Stack_Done;
            free(pt);
            free(str);
        }

END_OF_C_CODE

        # Place to store compiled inline code, create if non-existent
        $bin_store_dir = File::Spec->catfile( dirname(__FILE__), '.Inline' );
        mkdir $bin_store_dir unless -d $bin_store_dir;

        # Configure Inline to make taint mode work. Generates ugly warnings.
        if (${^TAINT}) {                # only when taint mode active
            require Inline;
            Inline->import(qw(untaint unsafe));
        }
    }

    # Compile code and store in the created directory (NOT in BEGIN block!)
    use Inline C => $c_code, directory => $bin_store_dir;
}

# Returns the universal gas constant in kcal/mol as defined in ViennaRNA.
sub get_gas_const {
    return $GASCONST;
}


# Convert pairtable and energy value to a formatted string. Pass structure and
# energy of structure. Excepts either a reference to a pairtable or a
# dot-bracket string as first arg.
sub printStrNrg
{
    my ($ptref, $nrg) = @_;
    # 0+nearest(0.01,$nrg) removes "negative zeroes", i.e. -0.00, because a
    # negative zero can only be created by negating a FLOAT zero, but
    # the result of  0+(-0.00) is an integer and therefore a "positive" zero
    # http://rosettacode.org/wiki/Extreme_floating_point_values#Perl
    return sprintf( "%s %6.2f", (ref($ptref) ? pt2str($ptref) : $ptref),
                                0+nearest(0.01,$nrg)
    );
}

# Converts a pair table to a dot-bracket string.
# Arguments:
#   pair_table: list or array ref to pairing table.
# Returns dot-bracket string of the structure.
sub pt2str( + ) { struct_string(@_) }
# # Consider using pt2str_c
# sub pt2str( + )
# {
#     my $ptref = shift;
#     my $form = '';
#     my $i = 0;
#     for my $j (@$ptref)
#     {
#         if( $j<0) { $form .= '.'; }
#         elsif( $i<$j) { $form .= '('; }
#         else { $form .= ')'; }
#         $i++;
#     }
#     return $form;
# }


# Generate pair table from dot-bracket string.
# Indexing starts at 0, -1 means unpaired (unlike ViennaRNA!). Returns an
# arrayref in scalar context, or a list of base pairs in list context.
sub str2pt
{
    my @str = split //, $_[0];
    my @pt;        # Pair table
    my @is;        # Index stack
    my $i=0;

    foreach my $c (@str) {
        for ($c) {
            when( '(') { push @is, $i; }
            when( ')') {
                croak "Closing paranthesis at position $i is unmatched."
                    if $#is==-1;
                my $j = pop @is;
                $pt[$i]=$j;
                $pt[$j]=$i;
            }
            when( '.') { $pt[$i]=-1; }
            default {
                croak "Structure contains invalid symbol '$c' at position $i "
                        . "(allowed are '(', ')' and '.').";
            }
        }
        $i++;
    }
    croak "Opening paranthesis at position $is[-1] is unmatched." if $#is!=-1;
    return @pt if wantarray;
    return \@pt;
}

# Normalize sequence, i.e. convert it to upper case, and delete any
# white-space characters.
# Arguments:
#   Sequence: sequence to normalize.
# Returns the normalized sequence.
sub normalize_seq ($) {
    my ($sequence) = @_;

    my $normalized_sequence = uc($sequence) =~ s{\s+}{}gr; # =~ s{T}{U}gr;

    return $sequence;
}

# Given a sequence string, structure string, or a pair table ref, return the
# sequence length.
sub seq_len ($) {
    my ($seq_or_ptref) = @_;
    my $seq_len
        = reftype $seq_or_ptref ? @$seq_or_ptref : length($seq_or_ptref);
    return $seq_len;
}

# Return true iff a valid sequence is given, that is one consisting solely of
# A, U, T, G, and C, either in upper or lower case. The empty sequence is not
# valid.
sub is_valid_seq ($) {
    my ($sequence) = @_;

    my $is_valid_seq = $sequence =~ m{ ^ [AUTGC]+ $ }ix;
    return $is_valid_seq;
}

# Check base pair validity.
# Pass sequence string and pair table reference. WARNING: This is slow and
# should be improved, cf. code.
# Returns 1 iff structure is valid.
sub bp_valid ( $+ ) { bpValid(@_) }             # snake case wrapper
sub bpValid ( $+ )
{
    #die "bpValid: Pass sequence string and pair table" unless @_==2;
    my @seq = split //, $_[0]; # convert to upper case
    my $ptref = $_[1];

    for( my $i=0; $i<@$ptref; $i++)
    {
        my $j = $ptref->[$i];
        if ($i < $j)    # $i is paired and pair has not been checked yet
        {
            # TODO a hash should be used for a constant-time comparison.
            unless( any {"$seq[$i]$seq[$j]" eq $_} @VALID_BP)
            {
                say STDERR "ERROR: bpValid: Invalid base pair ",
                           "'$seq[$i]$seq[$j]' at position ($i,$j).";
                return 0;
            }
        }
    }
    return 1;
}


# Parallely make passed (SORTED!) lists of equal length unique w.r.t. first list.
# More precisely, in the first list, every item is compared to the previous
# one, and if both are equal, it is discarded. If position i is discarded in
# the first list, it is also discarded in all other lists. The passed lists
# are modified in place.
sub parUniq
{
    my $fst = shift;
    #die "parUniq: First list has <2 elements or is undefined" unless @$fst>1;
    croak "parUniq: First list is undefined" unless defined($fst);
    for my $lst (@_) {
        croak "parUniq: Passed lists differ in their length" unless @$lst==@$fst;
    }
    return if @$fst<2;
    my $j = 1;                  # points to the last element not written to
    for my $i (1..$#$fst)
    {
        next if( $fst->[$i] ~~ $fst->[$i-1]);
        $fst->[$j] = $fst->[$i];
        for my $lst (@_) { $lst->[$j] = $lst->[$i]; }
        $j++;
    }

    # Discard superflous trailing elements
    splice @$_, $j foreach ($fst, @_);
}


# Parallely sort lists of equal length w.r.t. first list.
# More precisely, the first list is sorted numerically, and to the other
# lists, the same ordering is applied.
# The lists are sorted in place.
sub parSort
{
    my @tmp;
    my $fst = shift;
    croak "parSort: First list is undefined" unless defined($fst);
    return () unless (@$fst);     # First list contains no elements
    my @isort = sort { $fst->[$a] <=> $fst->[$b] } 0..$#$fst;
    foreach my $i (0..$#$fst) { push @tmp, $fst->[$isort[$i]]; };
    @$fst = @tmp;
    foreach my $lst (@_)
    {
        croak "parSort: Lists have different sizes" unless $#$fst == $#$lst;
        @tmp = ();
        foreach my $i (0..$#$fst) { push @tmp, $lst->[$isort[$i]]; };
        @$lst = @tmp;
    }
}


# Comparison function for two structure pair tables. Let struct1, struct2 be
# the two input structures. We define struct1 < struct2 iff
# sequence_length(struct1) < sequence_length(struct2) or
# sequence_length(struct1) == sequence_length(struct2) and
#     struct1 cmp struct2 < 0 (i.e. struct1 is lexicographically smaller)
# and use the total order induced by this relation.
# Arguments:
#   pt1: Pair table of first structure.
#   pt2: Pair table of second structure.
# Returns a comparison value, i.e. a negative value if pt1 < pt2, a positive
# value if pt1 > pt2, and zero if pt1 == pt2.
sub struct_cmp (++) {
    my ($pt1, $pt2) = @_;

    # If structures differ in length, the shorter one is smaller
    my ($last_index_pt1, $last_index_pt2) = ($#$pt1, $#$pt2);
    my $cmp_result = $last_index_pt1 <=> $last_index_pt2;
    return $cmp_result if $cmp_result;              # length not equal

    foreach my $i (0..$last_index_pt1) {
        # # Order: '.' < '(' < ')'
        # if ($pt1->[$i] == $pt2->[$i]) { $cmp_result= 0 }
        # elsif ($pt1->[$i] < 0)        { $cmp_result=-1 } # pt1=='.'
        # elsif ($pt2->[$i] < 0)        { $cmp_result= 1 } # pt2=='.'
        # elsif ($pt1->[$i] > $i)       { $cmp_result=-1 } # pt1=='(' & pt2==')'
        # elsif ($pt2->[$i] > $i)       { $cmp_result= 1 } # pt2=='(' & pt1==')'
        # else{ croak "oops, that's wrong"; }

        # Do lexicographic comparison: '(' < ')' < '.' which equals the
        # reversed numeric comparison of the two pair table entries
        $cmp_result = $pt2->[$i] <=> $pt1->[$i];

        return $cmp_result if $cmp_result;          # symbols differ
    }

    return 0;                                       # structures equal
}


# Check whether two structures are equal. For details, cf. struct_cmp.
# Arguments:
#   pt1: (Ref to) pair table of first structure.
#   pt2: (Ref to) pair table of second structure.
# Returns true value if the structures are equal, and a false value otherwise.
sub struct_eq (++) {
    return struct_cmp($_[0], $_[1]) == 0;
}


# Check whether a structure is smaller than another one. For details, cf.
# struct_cmp.
# Arguments:
#   pt1: (Ref to) pair table of first structure.
#   pt2: (Ref to) pair table of second structure.
# Returns true value if the first structure is smaller than the second one,
# and a false value otherwise.
sub struct_lt (++) {
    my ($pt1, $pt2) = @_;
    return struct_cmp($_[0], $_[1]) < 0;
}


# Compute the base pair distance of two structures given either as pair tables
# or dot-bracket strings.
sub bp_dist {
    my ($struct_a, $struct_b) = @_;

    # Convert dot-bracket strings to pair tables if required.
    my $pt_a = reftype $struct_a ? $struct_a : str2pt $struct_a;
    my $pt_b = reftype $struct_b ? $struct_b : str2pt $struct_b;

    confess 'Cannot compute base pair distance of structs with differnt lengths'
        unless @$pt_a == @$pt_b;

    my $bp_dist = 0;
    while (my ($i, $j1) = each @$pt_a) {
        my $j2 = $pt_b->[$i];
        next if $j1 == $j2;     # both unpaired or same bp
        $bp_dist++ if $j1 > $i; # opening bp in struct1 not present in struct2
        $bp_dist++ if $j2 > $i; # opening bp in struct2 not present in struct1
    }

    return $bp_dist;
}



# Extract a value from an option hash. If the option does not exist in the
# passed hash, a default value is returned if it was passed, otherwise undef.
# Arguments:
#   opt_ref:    Reference to a hash containing options.
#   opt_name:   String of the options name.
#   default:    Optional, default value if opt_name is not defined.
# Returns the value of the option if it exists in the opt_ref hash, otherwise
# a default value if it was passed, otherwise undef.
sub get_opt {
    my ($opt_ref, $opt_name, $default) = @_;
    $default = defined $default ? $default : undef;
    return exists $opt_ref->{$opt_name} ? $opt_ref->{$opt_name} : $default;
}


# Compute the partition function of each basin, normalized with the mfe.
# Arguments:
# sbmap_filehandle: Handle to a file that lists on each line a  structure,
#   its energy and its basin number; separated by whitespace. The structure
#   is not used and can be replaced by a non-whitespace dummy.
# temperature: Folding temperature in deg Celsius. [37]
# mfe: Optional, minimum free energy of sequence as a scalar. If undefined,
#   it is determined by reading the first line from the sbmap file (assumes
#   ordering by energy!). If a reference is provided, the mfe is determined
#   from sbmap_filehandle and the referenced scalar is updated with the
#   result. Otherwise, the given mfe is used.
# Returns a list of scaled basin partition functions, where the index
# corresponds to the basin number.
sub get_basin_partition_funcs {
    my ($sbmap_filehandle, $opt) = @_;
    my $temperature = get_opt($opt, 'temperature', 37);
    my $mfe_ref = get_opt($opt, 'mfe');

    if (not defined $mfe_ref or ref $mfe_ref) {     # find mfe
        # simply use energy value from first line...
        my $first_line = <$sbmap_filehandle>;
        my @fields     = split /\s+/, $first_line;
        $$mfe_ref      = $fields[1];

        seek $sbmap_filehandle, 0, 0;               # set fh position to start
    }
    # Make mfe a ref to the value it contains if it is a scalar
    $mfe_ref = do { my $tmp = $mfe_ref; \$tmp } unless ref $mfe_ref;

    #Second pass: compute partition function for each basin
    my @Z;    # partition function of basin i
    # BETA=1/RT. R = Boltzmann's constant, T = temperature in Kelvin
    my $BETA = 1/(get_gas_const()*celsius2kelvin($temperature));
    while( <$sbmap_filehandle>) {
        my @fields = split;
        my ($energy, $basin) = @fields[ 1..2 ];

        $Z[ $basin ] += exp $BETA*($$mfe_ref-$energy); # normalize w.r.t. mfe
    }
    return @Z;
}


# Given an RNA sequence of length n, compute an n x n probability matrix bpp
# storing the probabilities that nucleotide i is paired with nucleotide j in
# entry bpp[i,j]. This matrix is symmetric.
#
# Instead of a sequence, a ViennaRNA fold compound may be passed to allow
# adding constraints, etc. to the sequence.
#
# Arguments:
#   seq_or_fc:  Nucleotide sequence or ViennaRNA fold compound containing the
#               sequence.
#
# Returns reference to 2-dim array containing probabilities that index i pairs
# with index j.
sub base_pair_prob {
    my ($seq_or_fc) = @_;

    # Unless a fold compound has been passed, create a new one from the seq.
    my $fc = (blessed($seq_or_fc) || q{}) eq 'RNA::fold_compound'
             ? $seq_or_fc
             : RNA::fold_compound->new($seq_or_fc);

    # Required to make ViennaRNA compute partition function matrices
    $fc->pf();
    my $bpp = $fc->bpp();           # get raw base pair probability matrix

    # Remove empty first row / column
    shift @$bpp;
    shift @$_ foreach @$bpp;

    # Add symmetric entries to upper triangle matrix
    foreach my $i (0..($#$bpp-1)) {
        foreach my $j (($i+1)..$#$bpp) {
            $bpp->[$j][$i] = $bpp->[$i][$j];
        }
    }

    return $bpp;
}

# Convert a temperature value from degree Celsius to Kelvin.
sub celsius2kelvin ($) {
    return $_[0] + 273.15;
}


# Compute the Boltzmann weight Z of a given energy value E, i.e.
#       Z = exp( - E / (R * T) ),
# where T is the temperature and R is the universal gas constant.
# Arguments:
#   energy: energy of which to compute the Boltzmann weight, given in kcal/mol
# Optional arguments:
#   temperature: temperature at which to compute the Boltzman weight, given in
#       degree Celsius
# Returns the Boltzmann weight / factor.
sub nrg2boltz {
    my ($energy, $temperature) = @_;
    # Default temperature is 37 deg C
    $temperature = celsius2kelvin($temperature // 37);
    state $gas_const = get_gas_const();

    my $boltzmann_weight = exp( - $energy / ($gas_const*$temperature) );
    return $boltzmann_weight;
}


# Compute the energy value E associated with the given Boltzmann weight or
# partition function Z, i.e.
#       E = - R * T * ln Z
# where T is the temperature and R is the universal gas constant.
# Arguments:
#   boltzmann_weight: Boltzmann weight or partition function for which to
#       compute the associated energy value
# Optional arguments:
#   temperature: temperature at which to compute the Boltzman weight, given in
#       degree Celsius
# Returns the energy value.
sub boltz2nrg {
    my ($boltzmann_weight, $temperature) = @_;
    # Default temperature is 37 deg C
    $temperature = celsius2kelvin(defined $temperature ? $temperature : 37);
    state $gas_const = get_gas_const();

    return 'Inf' unless $boltzmann_weight;      # zero part func = inf energy

    my $energy = -$gas_const * $temperature * log $boltzmann_weight;
    return $energy;
}


# Returns an iterator that sequentially iterates through the passed array.
sub array_iter_factory {
    my $array_ref = shift;
    my $index = 0;

    return sub {
        if ($index < @$array_ref) {
            return $array_ref->[$index++];
        }
        else {
            return;
        }
    };
}


# Internal helper function used to verify integrity of optional arguments
# passed to a function. It simpy tests whether all passed options are found in
# the passed list of valid / allowed options and croaks when they're not.
# Args:
#   opt_ref: hash ref of options / optional arguments and their associated
#            value
#   valid_options: list of all allowed options
# Void. Croaks if non-allowed option is found.
sub _check_opts_valid {
    my ($opt_ref, @valid_options) = @_;

    # Create look-up table for valid options
    my %opt_is_valid = map {$_ => 1} @valid_options;

    # Verify that all options are valid, die if not
    foreach my $opt_name (keys %$opt_ref) {
        my $calling_func = (caller 1)[3];       # who called the check?
        croak "$calling_func(): invalid option '$opt_name'"
            unless $opt_is_valid{$opt_name};
    }
}


# Create a model details object which can be passed to ViennaRNA's fold
# compound constructor to set model details for the energy computation.
# Arguments: hash ref containing:
#   temperature: the RNA folding temperature in deg C [37]
#   dangles: dangling end model, either 0 (= no dangles), 1 , 2
#   (full/overdangle), or 3 (+ helix stacking) [2]
#   noLP: prohibit lonely base pairs (1 - only canonical structs, 0 - all) [0]
#         DOES NOT WORK AS EXPECTED for partition function calculations!
#         E.g. subopts will still contain some structures with lonely pairs,
#         and partition function will also not be exact.
#   pf_smooth: smooth energy parameters for partition function calculations [1]
#              Important to have differentiable energies when varying
#              temperature. Disable to get exact partition functions when
#              comparing with individual structures' Boltzmann weight (e.g.
#              for coverage analyses, probabilities of structures etc.)
#   compute_bpp: Compute base pair probabilities when doing partition function
#                folding. [1] Disable improve performance.
# Returns the model details object.
sub make_model_details {
    my ($opt_ref) = @_;

    # Verify option validity
    my @valid_opts = qw(dangles temperature noLP pf_smooth compute_bpp);
    _check_opts_valid $opt_ref, @valid_opts;

    # Create model details object and set options
    my $model_details = RNA::md->new();
    foreach my $opt_name (grep {defined $opt_ref->{$_}} keys %$opt_ref) {
        $model_details->{$opt_name} = $opt_ref->{$opt_name};
    }

    return $model_details;
}


# Create a ViennaRNA fold compound object from a sequence. If model options
# are given as a hash ref, create a model detail object and integrate settings
# into the fold compound.
# Arguments:
#   seq: RNA sequence to create compound for.
# Optional args, given as hash ref: cf. make_model_details()
# Returns ViennaRNA fold compound of the given sequence and settings.
sub make_fold_compound {
    my ($seq, $opt_ref) = @_;

    my $fold_compound;
    if ($opt_ref) {
        # Create model details object and set options.
        my $model_details = make_model_details $opt_ref;

        $fold_compound = RNA::fold_compound->new($seq, $model_details);
    }
    else {              # no options given, skip model_detail construction
        $fold_compound = RNA::fold_compound->new($seq);
    }

    return $fold_compound;
}


# Compute the minimum free energy (mfe) of a given sequence, considering
# optional arguments like temperature and dangling end model.
# Arguments:
#   seq: RNA sequence to compute mfe for.
#   opt_ref: Hash ref containing options, see make_model_details().
# In scalar context, returns the mfe of the sequence. In list context, both
# the mfe structure and the mfe itself are returned.
sub rna_mfe {
    my ($seq, $opt_ref) = @_;

    my ($mfe_struct, $mfe);

    if (defined $opt_ref) {     # set options using model detail object
        # Create fold compound object and compute mfe.
        my $fc  = make_fold_compound $seq, $opt_ref;
        ($mfe_struct, $mfe) = $fc->mfe;
    }
    else {                      # use flat interface for default params
        ($mfe_struct, $mfe) = RNA::fold($seq);
    }

    return wantarray ? ($mfe_struct, $mfe) : $mfe;
}


# Compute the structure of minimum free energy for a given sequence,
# considering optional arguments like temperature and dangling end model.
# Arguments:
#   seq: RNA sequence to compute mfe structure for.
#   opt_ref: Hash ref containing options, see make_model_details().
# Returns the mfe structure as dot-bracket string.
sub rna_mfe_struct {
    my ($seq, $opt_ref) = @_;
    my ($mfe_struct, $mfe) = rna_mfe($seq, $opt_ref);

    return $mfe_struct;
}


# Compute the free energy of an RNA secondary structure for a given sequence,
# considering optional arguments like temperature and dangling end model.
# Arguments:
#   seq: RNA sequence to compute energy for
#   struct: secondary structure to evaluate
#   opt_ref: hash ref containing options, see make_model_details().
# Returns the energy of the given structure.
sub rna_eval_structure {
    my ($seq, $struct, $opt_ref) = @_;

    # Create fold compound object and compute mfe.
    my $fc     = make_fold_compound $seq, $opt_ref;
    my $energy = $fc->eval_structure( $struct );

    return $energy;
}


# Computes the "partition energy" for a list of free energies. This is the
# energy related to the partition function by the Boltzmann transformation.
# More precisely, for a set of structures X and its partition function Z given
# by Z = sum Boltz( E(x) ) over x in X, where E(x) is the energy of x, its
# partition energy E_p and is computed such that Z = Boltz( E_p ).
# Arguments:
#   seq: Sequence for the given list of structures.
#   energy_iter: iterator that returns, with each call, the next secondary
#       structure's energy.
# Optional arguments:
#   scale_energy: Energy offset which is substracted from all structures'
#       energies before their Boltzmann weight is calculated. This avoids
#       numerical issues. Default: mfe of sequence
#   opt_ref: hash ref containing additional folding energy options:
#       temperature: folding temperature for Boltzmann weight conversions
# Returns the partition energy of the given list of free energies.
sub partition_energy_nrgiter {
    my ($seq, $energy_iter, $scale_energy, $opt_ref) = @_;

    # Scale huge partition Boltmann weights to avoid numerical issues.
    unless (defined $scale_energy) {
        $scale_energy = rna_mfe $seq, $opt_ref;
    }

    # Compute and sum Boltmann weights of all the structures.
    # require Math::BigFloat;
    # Math::BigFloat->accuracy(25);           # significant decimal digits
    # my $partition_func = Math::BigFloat->new(0);    # sum of Boltzmann weights
    my $partition_func = 0;                     # sum of Boltzmann weights
    my $temperature = $opt_ref ? $opt_ref->{temperature} : undef;
    while ( defined (my $struct_energy = $energy_iter->()) ) {
        $partition_func += nrg2boltz(
                               $struct_energy - $scale_energy,
                               $temperature,
                           );
    };

    # Convert back to energy and undo scaling
    my $partition_energy = boltz2nrg(
                                $partition_func,
                                $temperature,
                           )
                           + $scale_energy;
    return $partition_energy;
}


# Given a sequence and a struct iterator, computes the partition energy of the
# set of structures returned by the iter. Wrapper for
# partition_energy_structiter() that accepts a structure iter instead of an
# energy iter.
# Arguments:
#   seq: Sequence for which to compute the structure energies.
#   struct_iter: Iterator that returns the next structure for which to add up
#       the energies.
# Optional arguments:
#   scale_energy: Energy offset which is substracted from all structures'
#       energies before their Boltzmann weight is calculated. This avoids
#       numerical issues. Default: mfe of sequence.
#   opt_ref: hash ref containing additional folding energy options:
#       temperature: folding temperature for Boltzmann weight conversions
# Returns the scaled partition energy of the given structures.
sub partition_energy_structiter {
    my ($seq, $struct_iter, $scale_energy, $opt_ref) = @_;

    # Number of digits of Boltzmann weight for a given energy:
    # -1 * E / ( get_gas_const*celsius2kelvin(37)*log(10) )  +  1

    my $fc = make_fold_compound($seq, $opt_ref);
    my $energy_iter = sub {          # returns next structure's energy
        my $next_struct = $struct_iter->() or return;
        return $fc->eval_structure( $next_struct );
    };

    my $partition_energy
        = partition_energy_nrgiter $seq, $energy_iter, $scale_energy, $opt_ref;
    return $partition_energy;
}

# Wrapper for the partition_energy_iter function that can be called using a
# a sequence and an array of structures.
# Arguments:
#   seq: Sequence for which to compute the structure energies.
# Optional arguments:
#   opt_ref: Model detail options, cf. make_model_details().
#   structures: List of secondary structure dot--bracket strings for which to
#       compute the free energy. If empty, computes energy of entire ensemble.
# Returns the total energy covered by the passed structures. If none were
# passed, return the full ensemble energy.
sub partition_energy {
    my $seq = shift;

    # If next arg is a hash ref, interpret it as model details options hash.
    my $opt_ref = (@_ and reftype $_[0] eq reftype {}) ? shift : undef;

    # If no structure is passed, compute partition function of full ensemble.
    unless (@_) {
        my $ensemble_energy = ensemble_energy($seq, $opt_ref);
        return $ensemble_energy;
    }

    # Compute energies of passed structures.
    my $fc = make_fold_compound($seq, $opt_ref);
    my @energies = map { $fc->eval_structure($_) } @_;

    # Scale huge or tiny partition function to avoid numerical issues. Use the
    # mfe energy of the structure set to scale the largest Boltz weight to 1.
    my $scale_energy = min @energies;

    my $energy_iter = array_iter_factory \@energies;    # turn list into iter
    my $partition_energy
        = partition_energy_nrgiter $seq, $energy_iter, $scale_energy, $opt_ref;
    return $partition_energy;
}


# Compute ensemble / partition free energy. It does the same as
# partition_energy($seq), but accepts RNA model detail options to adjust
# temperature or dangling model.
# NOTE: Use parameter pf_smooth => 0 when using the ensemble energy together
# with energies of single structures (e.g. when computing probabilities).
# Smoothing would lead to errors!
# Arguments:
#   seq:         sequence for which the ensemble energy should be computed
# Optional args passed in hash ref: cf. make_model_details().
# The optional arguments default to their ViennaRNA default values.
# Returns the (full) ensemble free energy.
sub ensemble_energy {
    my ($seq, $opt_ref) = @_;

    my $ensemble_energy;

    if (defined $opt_ref) {
        # Compute ensemble energy using an RNA fold compound object.
        my $fc = make_fold_compound($seq, $opt_ref);
        $ensemble_energy = ($fc->pf)[1];
    }
    else {
        # Use simplified ("flat") interface if no model details are changed.
        $ensemble_energy = (RNA::pf_fold($seq))[1]
    }

    return $ensemble_energy;
}

# Compute the partial ensemble energy of an enumerated energy band above the
# mfe structure.
# Arguments:
#   seq:        Sequence for which the ensemble energy should be computed.
#   e_thresh:   Relative enumeration threshold (w.r.t. mfe) that is supposed
#               to be enumerated.
# Optional args passed in hash ref: cf. make_model_details().
# The optional arguments default to their ViennaRNA default values.
# NOTE: The subopt noLP (no lonely pair) option does NOT WORK THE SAME as the
# implementations in this module. Specifically, some structures may still
# contain lonely pairs.
# Returns the (partial) ensemble free energy of the enumerated energy band.
sub subopt_energy {
    my ($seq, $e_thresh, $opt_ref) = @_;

    # Disable unnecessary computations. Init opt ref if undef.
    $opt_ref //= {};
    $opt_ref->{compute_bpp} = 0;

    # Create fold compound object with adjusted model details.
    my $fc = make_fold_compound($seq, $opt_ref);

    # Use mfe to scale partition function.
    my $mfe = ($fc->mfe)[1];

    # Arguments to subopt: range needs to be deka cal integer.
    my $e_thresh_dekacal = int($e_thresh * 100);

    # Define callback closure that sums up partition function and store data
    # here.
    my $part_func = 0.;
    my $e_thresh_abs = $mfe + $e_thresh;
    my $sum_up_subopt_nrg = sub {
        # my ($struct, $nrg, $data) = @_;
        my ($struct, $nrg) = @_;            # could also use data

        # Last call has undefined struct (but not energy!), skip. Also skip if
        # energy is above threshold. This may happen with dangles=1|3 since
        # partition folding supports only dangles=0|2 and, with the first
        # setting, enumerated structures are re-evaluated AFTER the
        # backtracking by the subopt command.  Stand-alone RNAsubopt does this
        # filtering step, too.
        return unless defined $struct and $nrg <= $e_thresh_abs;

        # $nrg = sprintf '%.2f', $nrg;        # round as RNAsubopt does
        # say STDERR join ' :: ', @_;         # debug -- print all structures

        # Scale energies with mfe (such that mfe has Boltzmann weight of 1) to
        # avoid numerical issues. If this is not sufficient, use better approx
        # of ensemble energy instead of mfe struct for scaling (cf. RNAlib
        # manual).
        $part_func += nrg2boltz $nrg - $mfe, $opt_ref->{temperature};
    };

    # Generate subopts and let callback function sum up partition function.
    $fc->subopt_cb($e_thresh_dekacal, $sum_up_subopt_nrg);

    # Unscale again.
    my $subopt_nrg = $mfe + boltz2nrg $part_func, $opt_ref->{temperature};

    return $subopt_nrg;
}

# Seed the random number generator (Mersenne Twister) with the passed value.
# This is ONLY required to reproduce results, not before normal use.
# When passing undef, the generator is auto-seeded with a random integer. The
# used seed is returned. This is FAR LESS random than the random generator's
# internal auto-seeding, so really don't use this unless you need to.
# Arugments:
#   seed: The value to be passed to Math::Random::MT::Auto->srand. Can be a
#       single integer, an array ref containing integers etc. Passing undef
#       will auto-seed the generator with a single integer.
# Returns the used seed.
sub rand_seed {
    my ($seed) = @_;
    $seed //= irand();
    srand($seed);
    return get_seed();
}

# Generate a random RNA sequence consisting of G, C, A and U. Optionally, the
# exact GC content of the sequence can be specified. If not specified, it is
# chosen randomly.
# Arguments:
#   length: ~ of the sequence to be generated.
# Optional arguments:
#   gc_content: Expected GC content of the generated sequence. 'Expected'
#       means that this value is used as a probability for sampling either a G
#       or a C, so the actual GC content of the generated sequences may vary
# Returns a random sequence with the requested GC content.
sub rand_seq {
    my ($length, $gc_content) = @_;
    unless (defined $gc_content) {              # do not prescribe GC content
        state $symbols = [ qw( A U G C ) ];
        my @random_indices = map { int rand scalar @$symbols} 1..$length;
        my $sequence = join q{}, @$symbols[ @random_indices ];
        return $sequence;
    }

    my $gc_count = round $length * $gc_content;
    my $ac_count = $length - $gc_count;

    my @nucleotides;
    foreach (1..$length) {
        # Adjust GC prob with number of already printed GCs. This ensures the
        # requested gc frequency is i) met exactly (up to rounding), and ii)
        # that this frequency is, on average also achieved for each position
        # of the generated sequences, especially there is no GC bias towards
        # beginning or end of the sequence.
        if (rand() < $gc_count / ($gc_count+$ac_count)) {
            push @nucleotides, (rand() < 0.5 ? 'G' : 'C');
            $gc_count--;
        }
        else {
            push @nucleotides, (rand() < 0.5 ? 'A' : 'U');
            $ac_count--;
        }
    }

    my $sequence = join q{}, @nucleotides;
    return $sequence;
}


# Print a list containing (pair table ref, energy) pairs,
# each entry prefixed with an optionally passed prefix.
sub printStrNrgLst
{
    my @lst = @{$_[0]};
    my $prefix = $_[1];
    $prefix = "" unless defined($prefix);
    while( my @nxt = splice @lst, 0, 2)    # splice two elem in each iteration
    {
        my ($ptref, $nrg) = @nxt;
        say $prefix, printStrNrg($ptref, $nrg);
    }
}


# Perform gradient walk
# Pass seq as string, reference to pair table, and optionally reference
# to array where all encountered structures will be stored
# Returns pair table ref and energy of reached structure.
sub gradWalk( $+;$$$$ )
{
    my ($seqstr, $ptref, $walk, $verb, $noLP, $nrg) = @_;
    #my @pt = @$ptref;
    # die "gradWalk: Pass sequence string and pair table reference!"
    #     unless @pt && $seqstr;
    #say "seqstr=\n$seqstr, pt2str=\n",pt2str(@pt);
    $nrg = RNA::energy_of_struct( $seqstr, pt2str($ptref)) unless defined $nrg;
    while (1)
    {
        say "Next on walk:      ", printStrNrg($ptref, $nrg) if $verb;
        my @nb = genNeighbors( $seqstr, $ptref, 0, $noLP, 1, $nrg); # With noLP and shift moves
        if(@nb==0 || $nb[1] == $nrg)
        {
            say "Minimum reached" if $verb;
               return ($ptref, $nrg);
        }
        $ptref = $nb[0];
        $nrg = $nb[1];

        if( $walk)    # print structure to output file
        {
            say "   gradWalk: storing structure in passed array" if $verb;
            push @$walk, printStrNrg( $ptref , $nrg);
        }
    }
}


# Auxiliary sub that takes structure table reference and converts
# the structure pair table into the structure dot bracket string. This reduces
# the line noise that is caused by statements like
#   $seen{str2pt_c(@{$_->[0]})} = 1 foreach @struct
# to a more readable
#   $seen{ struct_string $_->[0] } = 1 foreach @struct
#
sub struct_string (+) {
    return pt2str_c( @{ $_[0] } );
}


# Perform single gradient walk step, i.e. a step to the neighbored structure
# of the lowest energy. To handle degenerate landscapes, if more than one
# structure has minimal energy, choose the smallest structure among them.
# Arguments:
#   sequence: sequence string
#   start_struct:  reference to pairtable of start structure
#   start_energy:  energy of start structure (optional, is computed if not
#                  defined)
#   opt:    option hash ref:
#       verb: flag, be verbose? [0]
#       noLP: flag, use moveset avoiding lonely pairs? [0]
#       shift: flag, use moveset using shift moves? [1]
# Returns the next structure of the grad walk and its energy or undef if the
# current structure is a local minimum.
sub grad_step ($+;$$) {
    my ($sequence, $current_struct, $current_energy, $opt) = @_;

    $current_energy = RNA::energy_of_struct($sequence,pt2str($current_struct))
        unless defined $current_energy;
    my $do_shift = defined $opt->{shift} ? $opt->{shift} : 1;

    my $energy_delta = 0;       # only return structures as good as input
    my $do_sort      = 0;       # do not sort returned neighbors
    my @neighbors
        = genNeighbors_pairs( $sequence, $current_struct, $energy_delta,
                              $opt->{noLP} || 0, $do_shift, $current_energy,
                              $do_sort
                            );

    # Leave structure undefined to recognize minima where no step is done
    my $min_struct;                 # = $current_struct;
    my $min_energy = $current_energy;

    foreach my $neighbor (@neighbors) {
        my ($neighbor_struct, $neighbor_energy) = @$neighbor;
        if ($neighbor_energy < $min_energy) {
            $min_energy = $neighbor_energy;
            $min_struct = $neighbor_struct;
        }
        elsif ($neighbor_energy == $current_energy) {   # saddle or plateau,
            # MUST come before next 'else if' because of initialization
            $min_struct = $neighbor_struct
                if $current_energy == $min_energy
                   and struct_lt $neighbor_struct, $current_struct;
        }
        elsif ($neighbor_energy == $min_energy) {
            # Keep lexicographically smallest
            $min_struct = $neighbor_struct
                if struct_lt $neighbor_struct, $min_struct;
        }
    }

    return defined $min_struct ? ($min_struct, $min_energy) : ();
}


# Perform a gradient walk from current structure down to local minimum of the
# current basin. Correctly handle degenerate landscapes, cf. grad_step
# TODO implement stop structure hash that allow shortcutting when doing many
# gradient walks from neighbored structures (e.g. in floodBasin)
# Arguments:
#   sequence: sequence string
#   start_struct:  reference to pairtable of start structure
#   start_energy:  energy of start structure (optional, is computed if not
#                  defined)
#   opt:    option hash ref:
#       walk_structs: array ref used to store structs (pt's) and energies of
#           structures encountered on walk between first and last struct if
#           they are required
#       verb: flag, be verbose? [0]
#       noLP: flag, use moveset avoiding lonely pairs? [0]
#       shift: flag, use moveset using shift moves? [1]
# Returns the last structure (pairtable) of the walk and its energy.
sub grad_walk($+;$$)
{
    my ($sequence, $start_struct, $start_energy, $opt) = @_;

    $start_energy = RNA::energy_of_struct($sequence, pt2str($start_struct))
        unless defined $start_energy;

    # Store encountered structures here
    my $walk_structs = $opt->{walk_structs} || [];

    my ($current_struct, $current_energy) = ($start_struct, $start_energy);
    my ($next_struct, $next_energy)
            = grad_step $sequence, $current_struct, $current_energy, $opt
        or return ($start_struct, $start_energy);   # already a min

    # Walk downward as long as energy reduces.
    while (defined $next_struct) {
        # Store structures of gradient walk
        ($current_struct, $current_energy) = ($next_struct, $next_energy);
        push @$walk_structs, ($current_struct, $current_energy);
        ($next_struct, $next_energy)
            = grad_step $sequence, $current_struct, $current_energy, $opt;
    }

    # Remove the last structure from the walk list and return it
    return (splice @$walk_structs, -2);
}

# Cycle (basin) exploration, algorithm description in summary_riboswitches!
# Flood cycle up to minh above starting structure's minimum and report all
# surrounding local minima reachable from transient structures found
# The following optional params are mainly used by waterfall:
# qedref:   Optional, pass ref to global hash of visited structures to re-use
#           this information across multiple calls.
# minNmin:  Optional, used with qedref. Only update qedref hash if at least
#           <minNmin> minima were found (allows re-visiting structures after
#           another minh has been passed.
# maxFlood: Optional, used with qedref and minNmin. If minh >= maxFlood,
#           update qedref in any case, even if not enough minima were found
# ofh:      Output file handle, write structures and their energies found in
#           this cycle to a file.
# noLP      Use moveset avoiding lonely pairs by removing/adding lonely
#           stacks of size 2 if a single pair indel would be lonely. [1]
sub floodCycle( $$$;$$$$$)
{
    my $verb = 0;
    my ($seqstr, $startptref, $minh, $qedref, $minNmin, $maxFlood, $ofh, $noLP) = @_;
    $qedref = {} unless defined $qedref;     #Reference to hash of structures already queued
    $minNmin = -1 unless defined $minNmin;
    $maxFlood = -1 unless defined $maxFlood;
    $noLP = 1 unless defined $noLP;
    my $shift = 1;
    my %qedloc;         # Local hash for queued structures
    my %minloc;         # Local hash for found minima
    my (@min, @minNrg);    #list of surrounding minima
    my @write;          # List of strings to write to ofh
    my @write_maybe;    # List of strings to write to ofh if not yet written
    assert $startptref && $seqstr;

    my $startptrefstr = pt2str(@$startptref);
    my $startnrg = RNA::energy_of_struct( $seqstr, $startptrefstr);
    my @Q = ($startptref, $startnrg); # Structure queue
    $qedloc{ $startptrefstr } = 1;
    $minloc{ $startptrefstr } = 1;

    my $maxnrg = $startnrg + $minh;
    my $nstr = 0;

    # Generate list of transient states that lead to lower minima than start structure.
    while( @Q)
    {
        say "\tExploring structure ", ++$nstr, "... |Q|=", @Q/2 if $verb;
        my $nrg   = pop @Q;
        my $ptref = pop @Q;
        say printStrNrg( $ptref, $nrg ) if $verb;
        assert defined $nrg;

        my @walk;        # store structures encountered on walk here
        my ($gradptref, $gradnrg)
            = gradWalk( $seqstr, $ptref, defined($ofh) ? \@walk : undef,
                        0, $noLP, $nrg);
        if($gradnrg<$startnrg)    # landed below start structure, found transient structure!
        {
            my $gradstr = pt2str @$gradptref;
            next if $minloc{ $gradstr };        # already found this minimum
            $minloc{ $gradstr } = 1;

            push @write_maybe, map {split} @walk;     # store encountered structures
            say "\t    Found trans struct" if $verb;
            push @min, $gradptref;
            push @minNrg, $gradnrg;
        }
        else
        {
            my @N = genNeighbors( $seqstr, $ptref, $maxnrg-$nrg, $noLP, $shift, $nrg);
            while( my ($nxt_ptref, $nxt_nrg) = splice @N, 0, 2)    # splice two elem in each iteration
            {
                my $nstrstring = pt2str(@$nxt_ptref);

                # Do not explore known structures again
                next if $qedloc{ $nstrstring } || $qedref->{ $nstrstring };
                $qedloc{ $nstrstring } = 1;    # mark as queued

                # Keep structure and nrg to write it to file later in last iteration
                say("\t    Pushing something") if $verb;
                push @write, $nstrstring, $nxt_nrg if defined $ofh;
                push @Q, $nxt_ptref, $nxt_nrg;
            }
        }
    }

    # Sort and de-dupe return list
    assert(@min == @minNrg);
    my $nmin = @min;
    parSort( \@minNrg,\@min);    # Return a sorted list of minima
    parUniq( \@min,\@minNrg);    # Remove duplicate items
    say "   Removed ", $nmin-@min, " duplicate minima" if $nmin>@min && $verb;
    say "   Update? |min|=", scalar @min, "  maxFlood=$maxFlood  minh=$minh " if $verb;
    if( $minNmin>=0 && @min>=$minNmin || $maxFlood>=0 && $maxFlood-$minh<0.001)
    {    # found enough minima or reached max flood lvl
        say "   Updating global qed hash" if $verb;
        @{$qedref}{ keys %qedloc } = values %qedloc;
        if( defined $ofh)
        {
            say "   floodCycle: Writing structures to output file" if $verb;
            while (my ($strstring, $nrg) = splice @write_maybe, 0, 2) {
                # Don't write structures twice
                next if $qedref->{ $strstring };

                say {$ofh} printStrNrg($strstring, $nrg);
            }

            while (my ($strstring, $nrg) = splice @write, 0, 2) {
                say {$ofh} printStrNrg($strstring, $nrg);
            }
        }
    }
    return () unless @min;
    return map { ($min[$_], $minNrg[$_]) } 0..$#min;
    # my @ret;
    # for( my $i=0; $i<@min; $i++)
    #     { push @ret, $min[$i], $minNrg[$i]; }
    # return @ret;
}


# Get the best (i.e. energetically lowest) direct path connecting
# the start structure with the target structure. The search is heuristic
# and keeps at most n structures at a time. The search space can be
# restricted to canonical paths (noLP, no lonely pairs). Shift moves
# are also supported.
# Arguments:
#   seq: Sequence string.
#   start: Start structure as dot-bracket string.
#   target: Target structure as dot-bracket string.
#   n: Keept at most n structures during search.
#   noLP: Do not allow lonely pairs, i.e. search canonical paths [0]
#   shift: Allow shift moves.  [0]
# Returns a list of structure strings, including their energies, denoting the
# found optimal path.
sub bestDirectPath( $$$;$$$)
{
    my ( $seq, $start, $target, $n, $noLP, $shift,) = @_;
    $noLP  //= 0;
    $shift //= 0;
    my ( $pta, $ptb) = ( scalar str2pt($start), scalar str2pt($target));
    my ( $is_val_pta, $is_val_ptb)
        = (bpValid( $seq, $pta), bpValid( $seq, $ptb));
    croak "Invalid basepair found in "
        . (!$is_val_pta ? "start" : "target")
        . " structure."
        unless $is_val_pta and $is_val_ptb;

    my $nCur = 1;            # Iteratively increase structures kept to n
    my $maxNrg   = undef;
    my $bestPath = undef;
    my $opt = { noLP=>$noLP, shift=>$shift, nreturn=>1, noDupes=>1 };

    while( 1) {
        $opt->{maxNrg} = $maxNrg;        # Update current upper energy bound
        my $newPath = (directPaths( $seq, $pta, $ptb, $nCur, $opt))[0];
        if( $newPath) {    # Found a better path
            $bestPath = $newPath
                if !defined($bestPath)
                    || $bestPath->{saddle}>$newPath->{saddle};
            $maxNrg = $bestPath->{saddle};
        }
        last if $nCur==$n;    # did search with $n paths
        $nCur = min(2*$nCur, $n);
    }
    croak 'Oops, no path found, did you request noLP on a non-canonical struct?'
        unless $bestPath;                # none found
    return path2str( $bestPath );
}


# Get the best (i.e. energetically lowest) direct path connecting
# the start structure with the target structure and return its saddle
# energy. The search is heuristic and keeps at most n structures at a
# time. The search space can be restricted to canonical paths
# (noLP, no lonely pairs). Shift moves are also supported.
# Arguments:
#   seq: Sequence string.
#   start: Start structure as dot-bracket string.
#   target: Target structure as dot-bracket string.
#   n: Keept at most n structures during search.
#   noLP: Do not allow lonely pairs, i.e. search canonical paths.
#   shift: Allow shift moves.
# Returns the saddle (i.e. highest) energy of the best (i.e. energetically
# lowest) path found, rounded to two digits.
sub bestDirectSaddle( $$$;$$$)
{
    my ( $seq, $start, $target, $n, $noLP, $shift,) = @_;
    my ( $pta, $ptb) = ( scalar str2pt($start), scalar str2pt($target));
    my ( $vala, $valb) = (bpValid( $seq, $pta), bpValid( $seq, $ptb));
    croak "Invalid basepair found in " . (!$vala ? "start" : "target") . " structure."
        unless $vala and $valb;
    my $saddle = ${ (directPaths( $seq, $pta, $ptb, $n,
        { noLP=>$noLP, shift=>$shift, nreturn=>1, noDupes=>0 }))[0] }{saddle};
    return nearest( 0.01, $saddle); # round saddle energy to two digits
}


# Compute direct paths with low saddles from structure pta to ptb (pass pairtables).
# Keep at most n entries during search.
# Arguments:
#   seqstr: Sequence as string.
#   pta: Pair table of start structure.
#   ptb: Pair table of target structure.
#   n: Keep at most n structures during search.
#   opt: Optional hash for options, see below.
# Options:
#   shift: Allow shift moves. [1]
#   noLP: Only search canonical paths. [1]
#   verb: Be verbose.
#   nreturn: Return only the best <nreturn> paths. [Inf]
#   noDupes: No duplicates. [0]
#   maxNrg: Maximum allowed energy in the path. [Inf]
# Returns a list of the paths found.
sub directPaths( $++$;+ )
{
    my ($seqstr, $pta, $ptb, $n, $opt) = @_;
    my $shift   = defined( $opt->{shift})   ? $opt->{shift}   : 1;
    my $noLP    = defined( $opt->{noLP})    ? $opt->{noLP}    : 1;
    my $verb    = defined( $opt->{verb})    ? $opt->{verb}    : 1;
    my $nreturn = defined( $opt->{nreturn}) ? $opt->{nreturn} : -1; # return only best $nreturn
    my $noDupes = defined( $opt->{noDupes}) ? $opt->{noDupes} : 0;
    my $maxNrg  = defined( $opt->{maxNrg})  ? $opt->{maxNrg}  : undef;

    eval {
        require Heap::Priority;
    } or croak 'directPaths() requires the PATCHED version of Heap::Priority';

    # Build path of length 0 from $pta and add it to pathlist
    my @paths =  (genInitialPath( $seqstr, $pta, $ptb)); # reference to path lists
    my $changed = 1;

    while ($changed) {
        $changed = 0;
        @paths = continuePaths( \@paths, $ptb, \$changed, $seqstr,
                                { shift=>$shift, noLP=>$noLP, n=>$n,
                                  noDupes=>$noDupes, maxNrg=>$maxNrg
                                }
                              );    # Add new paths

        assert( [ sort { $a->{saddle} <=> $b->{saddle} } @paths ] ~~ @paths);
        assert @paths <= $n;
    }

    if ($nreturn > 0) {
        return @paths[0 .. min( $nreturn-1, $#paths)];
    }
    else {
        return @paths;
    }
}


# Elongate any path in the passed pathlist by any possible move
# paths: ref to list of current paths to elongate
# ptb: ref to pair table of target structure
# changed: ref to a scalar to signal that at least one path was
#          elongated
# seqstr: sequence string
# opt: ref to hash with optional parameters
#     shift: allow shift moves
#     noLP: generate only canonical structures / paths
#     n: only return the n best neighbors
#     noDupes: flag, set to prevent returning structures with equal
#              pair tables. Works only with n>0.
sub continuePaths( ++$$;+ )
{
    my ($paths, $ptb, $changed, $seqstr, $opt ) = @_;
    my $shift = defined( $opt->{shift}) ? $opt->{shift} : 1;
    my $noLP = defined( $opt->{noLP}) ? $opt->{noLP} : 1;
    my $verb = defined( $opt->{verb}) ? $opt->{verb} : 1;
    my $n = defined( $opt->{n}) ? $opt->{n} : -1;
    my $noDupes = defined( $opt->{noDupes}) ? $opt->{noDupes} : 0;
    my (@newPaths, %cued);     # store new paths, remember queued structures
    my $maxNrg = defined( $opt->{maxNrg}) ? $opt->{maxNrg} : '-Inf';
    my $boundNrg = defined( $opt->{maxNrg});    # initial upper nrg bound passed?
    my $cue;
    $cue = new Heap::Priority, $cue->raise_error(2) if $n>0; # create new queue with error reporting
    foreach my $path (@$paths) { # Continue each path on list
        my @neigh;
        if( $path->{dist2target}>0 ) { # Not yet at target
            $$changed = 1; # We're generating new neighbors

            @neigh = genDirectedNeighbors( $path, $ptb, $seqstr,
                                             {
                                                 shift      =>  $shift,
                                                 noLP       =>  $noLP,
                                                 noDupes    =>  $noDupes,
                                                 maxNrg     =>
                                                   ($boundNrg or $n==$cue->count)
                                                         ? $maxNrg
                                                         : undef,
                                             }
                                         );
        }
        else { # we are already at target structure
            @neigh = ($path);
        }

        if( $n>0) {     # keep only best n paths using a priority queue
            for my $newpath (@neigh) {
                my $sad = $newpath->{saddle};
                if( (my $firstn = $cue->count<$n) or $sad<$maxNrg) { # fill cue; then push when better
                    my $str = ''; $str = pt2str_c(@{$newpath->{last}}) if $noDupes;
                    if( $noDupes and exists $cued{$str}) { # structure currently cued
                        my $old = $cued{$str};
                        if( $sad < (my $oldsad=$old->{saddle})) { # found better path, replace!
                            $cue->delete_item($old, $oldsad); # some issue with floats
                            $cue->add( $newpath, $sad);
                            assert( scalar grep($newpath, $cue->get_list));
                            $cued{$str} = $newpath;    # update hash with new reference
                            $maxNrg = $cue->next_level; #update in case exchanged==next elem
                        }
                    }
                    else { # structure is no dupe or don't care
                        $cue->add( $newpath, $sad);
                        assert( scalar grep( $newpath, $cue->get_list));
                        if( $firstn) {    # add, count, hash, do not pop anything
                            $cued{$str} = $newpath if $noDupes;
                        }
                        elsif( $noDupes) {    # mark structure as cued, unmark popped element
                            my $popstr = pt2str_c(@{${$cue->pop}{last}});
                            #say STDERR "popstr: $popstr";
                            assert( exists $cued{$popstr});
                            delete $cued{$popstr};
                            $cued{$str} = $newpath;
                        }
                        else {
                            $cue->pop;              # drop last path, we found a better one
                        }
                        $maxNrg = $cue->next_level; # update current max nrg
                    }
                }
                assert( $cue->count == $cue->count);
            }
        }
        else {                                      # do not filter
            push @newPaths, @neigh;
        }
    }
    @newPaths = reverse $cue->get_list if $n>0;
    return @newPaths;
}


# Elongate the passed path by adding all valid target-directed neighbors.
# Returns a list of the elongated paths.
sub genDirectedNeighbors ( ++;+)
{
    my ($path, $ptb, $seqstr, $opt) = @_;
    my $shift = defined( $opt->{shift}) ? $opt->{shift} : 1;
    my $noLP = defined( $opt->{noLP}) ? $opt->{noLP} : 1;
    my $verb = defined( $opt->{verb}) ? $opt->{verb} : 1;
    my $maxNrgDif = defined( $opt->{maxNrg}) ? $opt->{maxNrg}-$path->{nrg} : undef;
    my $newNrgDif = undef;
    my @neigh;

    my $ind = -1;
    # say "last: \n01234567890123456789012345\n", pt2str $path->{last},
    #      " add: @{$path->{add}} rem: @{$path->{rem}}" ;
    for my $i (@{$path->{add}}) # add basepairs
    {
        $ind++;
        next if $i<0;    # already added bp
        my $j = $ptb->[$i];
        my $pt = $path->{last};
        if( $pt->[$i]<0 and $pt->[$j]<0 and isNonCrossing( $pt, $i, $j)) # i and j unpaired
        {
            if( $noLP and isLonely( $pt, $i, $j))
            {
                if( $ind<$#{$path->{add}} and $path->{add}[$ind+1]==$i+1
                        and $ptb->[$i+1]==$j-1 # consecutive pair is chi-direct
                        and $pt->[$i+1]<0 and $pt->[$j-1]<0 # inner pair unpaired
                        and isInsLonely($pt,$i+1,$j-1))
                {
                    if( defined $maxNrgDif) {    # pre-compute energy dif of double move
                        $newNrgDif = nrg_of_move($seqstr,$pt,$i,$j);
                        addBp( $pt, $i, $j);
                        $newNrgDif += nrg_of_move($seqstr,$pt,$i+1,$j-1);
                        remBp( $pt, $i, $j);
                    }
                    if( !defined($maxNrgDif) or $newNrgDif<$maxNrgDif)
                    {
                        my $newPath = clone $path;
                        updatePath( $newPath, $i, $j, $ind, $seqstr,
                            { isAdd=>1, isDoubleMove=>1, newNrgDif=>defined($newNrgDif) ? 0 : undef } );
                        updatePath( $newPath, $i+1, $j-1, $ind+1, $seqstr, { isAdd=>1, newNrgDif=>$newNrgDif });
                        push @neigh, $newPath;
                    }
                }

            }
            else # non-lonely or don't care, add as usual
            {
                if( !defined($maxNrgDif) or # don't perform move/clone if energy too high
                        ($newNrgDif=nrg_of_move($seqstr,$pt,$i,$j))<$maxNrgDif ) {
                    my $newPath = clone $path;
                    updatePath( $newPath, $i, $j, $ind, $seqstr, { isAdd=>1, newNrgDif=>$newNrgDif });
                    push @neigh, $newPath;
                }
            }
        }
    }

    $ind = -1;
    for my $i (@{$path->{rem}}) # rem basepairs
    {
        $ind++;
        next if $i<0;    # already removed bp
        my $pt = $path->{last};
        my $j = $pt->[$i];
        my ($ilg, $olg) = (insGrowLonely($pt,$i,$j), outGrowLonely($pt,$i,$j));
        next if $noLP and $olg; # Only handle inside-lonely-growing pairs
        if( $noLP and ($ilg or $olg))
        {
            if( isOutLonely($pt,$i,$j) and $ilg
                    and $ind<$#{$path->{rem}} and $path->{rem}[$ind+1]==$i+1)
            {
                if( defined $maxNrgDif) {    # pre-compute energy dif of double move
                    $newNrgDif = nrg_of_move($seqstr,$pt,$i,$j);
                    remBp( $pt, $i, $j);
                    $newNrgDif += nrg_of_move($seqstr,$pt,$i+1,$j-1);
                    addBp( $pt, $i, $j);
                }
                if( !defined($maxNrgDif) or $newNrgDif<$maxNrgDif) {
                    my $newPath = clone $path;
                    updatePath( $newPath, $i, $j, $ind, $seqstr,
                        { isAdd=>0, isDoubleMove=>1, newNrgDif=> defined($newNrgDif) ? 0 : undef });
                    updatePath( $newPath, $i+1, $j-1, $ind+1, $seqstr, { isAdd=>0, newNrgDif=>$newNrgDif });
                    push @neigh, $newPath;
                }
            }
        }
        else # not a lonely pair or don't care, remove
        {
            my $newPath = clone $path;
            if( !defined($maxNrgDif) or # don't perform move/clone if energy too high
                    ($newNrgDif=nrg_of_move($seqstr,$newPath->{last},$i,$j))<$maxNrgDif )
            {
                updatePath( $newPath, $i, $j, $ind, $seqstr, { isAdd=>0, newNrgDif=>$newNrgDif });
                push @neigh, $newPath;
            }
        }
    }

    if( $shift) # Perform shift moves
    {
        $ind = -1;
        for my $i (@{$path->{rem}})
        {
            $ind++;
            next if $i<0;    # already removed bp
            my $pt = $path->{last};
            my $j = $pt->[$i];
            my ($ib, $jb) = ($ptb->[$j], $ptb->[$i]);
            while ( $ib>=0 or $jb>=0 ) # shift i or shift j
            {
                my ($lft, $rgt);    # new bp that will be inserted
                if( $ib>=0) { # get order of indices of new basepair
                    ($lft, $rgt) = $ib<$j ? ($ib,$j) : ($j,$ib); $ib=-1; }
                else { # $jb>=0
                    ($lft, $rgt) = $i<$jb ? ($i,$jb) : ($jb,$i); $jb=-1; }
                if( $noLP) {    # avoid creation of lonely pairs if we care
                    next if growLonely($pt,$i,$j) or isLonely($pt,$lft,$rgt);
                }
                remBp( $pt, $i, $j);    # Remove bp for crossing check
                my $nonCross = isNonCrossing($pt,$lft,$rgt) && $pt->[$lft]<0 && $pt->[$rgt]<0;
                addBp( $pt, $i, $j);
                if( $nonCross )
                {
                    my $newPath = clone $path;
                    if( defined $maxNrgDif) {    # pre-compute energy dif of double move
                        $newNrgDif = nrg_of_move($seqstr,$pt,$i,$j);
                        remBp( $pt, $i, $j);
                        $newNrgDif += nrg_of_move($seqstr,$pt,$lft,$rgt);
                        addBp( $pt, $i, $j);
                    }
                    if( !defined($maxNrgDif) or $newNrgDif<$maxNrgDif)
                    {
                        updatePath( $newPath, $i, $j, $ind, $seqstr,
                            { isAdd=>0, isDoubleMove=>1, newNrgDif=> defined($newNrgDif) ? 0 : undef });
                        my $addInd = -1;
                        for my $add (@{$path->{add}}) { # get index of bp addition
                            $addInd++; last if $lft==$add; }
                        assert $addInd<@{$path->{add}};
                        updatePath( $newPath, $lft, $rgt, $addInd, $seqstr, { isAdd=>1, newNrgDif=>$newNrgDif });
                        push @neigh, $newPath;
                    }
                }
            }
        }
    }

    return @neigh;
}


# Add a new move to a path, i.e. update structure, nrg, saddle etc
# i,j: Basepair that is added or removed
# ind: index of i in the path's add / rem list
# seqstr: sequence
# isAdd: Flag, if set, add bp (i,j), else remove it
# isDoubleMove: Flag, if set, the current move is the first one of a double move
#               i.e. lonely stack, shift
# newNrg: if energy of the next move was already computed, use it.
sub updatePath( +$$$;+)
{
    my ( $newPath, $i, $j, $ind, $seqstr, $opt) = @_;
    my $isAdd = defined($opt->{isAdd}) ? $opt->{isAdd} : 1;
    my $isDoubleMove = defined($opt->{isDoubleMove}) ? $opt->{isDoubleMove} : 0;
    my $newNrgDif = defined($opt->{newNrgDif}) ? $opt->{newNrgDif} : undef;
    assert( ($newPath->{last}[$i]<0) == $isAdd);
    # say "01234567890123456789";
    #say pt2str $newPath->{last}, " i=$i j=$j isAdd=$isAdd";

    $newPath->{nrg} += defined($newNrgDif) ? $newNrgDif :
                             nrg_of_move( $seqstr, $newPath->{last}, $i, $j); # update energy
    push @{$newPath->{nrgs}}, $newPath->{nrg};
    if( $isAdd)
    {
        addBp( $newPath->{last}, $i, $j);
        $newPath->{add}[$ind] = -1;  # mark bp as added
    }
    else
    {
        remBp( $newPath->{last}, $i, $j);
        $newPath->{rem}[$ind] = -1;  # mark bp as added
    }
    push @{$newPath->{moves}}, [ $i, $j ];
    $newPath->{dist2target}--;
    # Don't set saddle height to nrgs of intermediate structures
    $newPath->{saddle} = max( $newPath->{saddle}, $newPath->{nrg}) unless $isDoubleMove;
    $newPath->{moveInd}++;    # increment move index to point to this move
    push @{$newPath->{doubleMoves}}, $newPath->{moveInd}+1 if $isDoubleMove;

    return $newPath;
}


sub genInitialPath( $$$)
{
    my ($seqstr, $pta, $ptb) = @_;
    my %path = (
        last => $pta,    # last structure pt of current path
        nrg => RNA::energy_of_struct( $seqstr, pt2str $pta), # nrg of current structure
        nrgs => [],     # store energies of all encountered structures => |nrgs| = |moves| + 1
        moves => [],
        moveInd => -1,    # Index of last added move, -1 if move list is empty
        # these moves together with preceding one are 2bp moves (these are moves list indices!)
        doubleMoves => [],
        add => [],
        rem => [],
        dist2target => -1,    # basepair distance
        saddle => 0, # highest nrg of path
    );
    $path{saddle}=$path{nrg};
    push @{$path{nrgs}}, $path{nrg};

    my $i = 0;
    for my $j( @$pta)    # Get all bp's that need to be added or removed
    {
        my $jb = $ptb->[$i];
        if( $j<0) # a unpaired
        {
            push @{$path{add}}, $i if $jb>$i; # ib paired and left end of bp
        }
        elsif( $i>$j) # a paired, closing bp
        {
            push @{$path{add}}, $i if $i<$jb;
        }
        else # a paired, opening bp
        {
            if( $jb<0) # b unpaired
            {
                push @{$path{rem}}, $i;
            }
            elsif( $j!=$jb) # different bp in target
            {
                push @{$path{rem}}, $i;
                push @{$path{add}}, $i if $i<$jb;
            }
        }
        $i++;
    }
    assert defined $path{add};
    assert defined $path{rem};
    $path{dist2target} = @{$path{add}} + @{$path{rem}};

    return \%path;
}


# Checks whether (i,j) causes a crossing when inserted into pt.
# i and j are assumed to be unpaired!
sub isNonCrossing( +$$)
{
    my ($pt, $i, $j) = @_;
    for( my $k=$i+1; $k < $j; $k++)
    {
        #say "NonCrossing: k=$k
        my $l = $pt->[$k];
        if( $k<$l) # skip substructure
        {
            $k=$l;
            return 0 if $j<$l; # j lies in substructure
        }
        elsif( $l>=0 ) # end of current loop
        {
            return 0;
        }
    }
    return 1;
}


# For the passed path, return a list of dot-bracket strings of visited structures
# from start to target.
sub path2str ( +)
{
    my $path = shift;
    # return if not $path;

    # Reconstruct structures from last pairtable in reversed order, then
    # reverse result again
    my ($str, $lastStr) = ( undef, $path->{last});      # current and last structure
    my @str = ($lastStr);                               # list of structures to return

    # Fix indices for reversed move list
    my @doubleMoves = map { $path->{moveInd}-$_ } reverse @{$path->{doubleMoves}};
    my @nrgsCpy = @{$path->{nrgs}};                     # copy energies to be non-destructive
    my @nrgs = ( pop @nrgsCpy);
    my $ind =0; # index of current move (list reversed!)
    # state variables for double move
    my ($nxtDoubleMove, $isSecondDouble) = (@doubleMoves ? shift @doubleMoves : -1, 0);
    for my $move ( reverse @{$path->{moves}} )
    {
        my ($i, $j) = @$move;
        $str = clone $lastStr unless $isSecondDouble;
        $isSecondDouble=0;          # reset state
        if( $str->[$i]<0) {
            addBp( $str, $i, $j);
        } else {
            remBp( $str, $i, $j);
        }

        if( $ind == $nxtDoubleMove) # double move, don't push yet
        {
            $nxtDoubleMove = @doubleMoves ? shift @doubleMoves : -1; # get next double move
            $isSecondDouble = 1;    # don't clone in next iteration
            pop @nrgsCpy;           # throw away energy of transient structure
        }
        else                        # no double move, push structure
        {
            push @nrgs, pop @nrgsCpy;
            push @str, $str;
        }

        $lastStr = $str;
        $ind++;
    }

    return reverse map { printStrNrg( $_, shift @nrgs) } @str;
}


######################################
##  Classical loneliness & moveset  ##
######################################


# Wrapper for RNA::energy_of_move from ViennaRNA to efficiently compute the
# energy difference of a single bp move.
# Arguments:
#   seqstr: sequence string
#   ptref:  reference to pairtable of struct
#   i, j:   base pair that WILL BE BUT IS NOT YET added/removed to/from
#           pairtable
# Returns energy difference that would arise when adding / removing this bp.
#
# Notes:
#   - ViennaRNA's pairtables are 1-based, in this package they are 0-based.
#   - ViennaRNA tries to ADD (i,j) if i,j>0, and to REMOVE (i,j) if i,j<0.
#   - ViennaRNA requires i<j or the package could die
sub nrg_of_move( $+$$ )
{
    my ($seqstr, $ptref, $i, $j) = @_;
    ($i, $j) = ($j, $i) if $j < $i;             # make (i,j) ordered
    $i++;                                       # make 1-based
    $j++;
    if ($ptref->[$i-1] >= 0) {                  # $i is paired: negate (i,j)
        $i *= -1;
        $j *= -1;
    }
    return RNA::energy_of_move( $seqstr, pt2str_c(@$ptref), $i, $j);

    # return RNA::energy_of_move( $seqstr, pt2str_c(@$ptref),
    #             $ptref->[$i]<0 ? $i+1 : -$i-1, $ptref->[$j]<0 ? $j+1 : -$j-1);
}


# Generate all neighbors and filter w.r.t. their energy difference to the
# current structure.
# Arguments:
#   seq: Sequence string.
#   pair table ref: Pair table of structure.
# Optional arguments:
#   minh:  Maximal energy difference of neighbors. For the energy E of the
#          current structure, only return structures with energy at most
#          E + minh. [Infinity]
#   noLP:  Skip neighbors containing lonely pairs? [1]
#   shift: Perform shift moves? [1]
#   nrg:   Energy of current structure. Re-computed unless provided.
#   sort:  Sort the returned structure list by energy. [1]
#   do_energies: Compute energy of neighbors? If not, all energy values are
#                set to 0 and options minh, nrg and sort are disabled. [1]
# Returns (sorted?) list of pairs (\@pt,$nrg) (ref to pair table, energy)
# resembling the neighbor structures fulfilling the energy criterion.
sub genNeighbors
{
    my $verb = 0;        # be verbose (for debugging)
    my $noLP = 1;        # avoid generation of structs containinglonely pairs
    my $shift = 1;       # perform shift moves
    my $sort = 1;        # sort the returned list?
    my $do_energies = 1; # should we compute move energies
    my $minh = 'Inf';    # maximum distance from start struct

    my $seqstr = $_[0];
    my @seq = split //, $seqstr;
    my @pt = @{$_[1]};
    croak "genNeighbors: Pass sequence string and pair table reference!"
        unless @pt && $pt[0] && $pt[1] && $seqstr;
    $minh        = $_[2] if defined $_[2];
    $noLP        = $_[3] if defined $_[3];
    $shift       = $_[4] if defined $_[4];
    my $nrg      = $_[5];
    $sort        = $_[6] if defined $_[6];
    $do_energies = $_[7] if defined $_[7];

    my $nrg_of_move_sub = \&nrg_of_move;
    unless ($do_energies) {             # disable energy computation
        $minh = 'Inf';                  # don't filter any structs
        $nrg  = 0;                      # don't compute mfe
        $sort = 0;                      # don't sort crap energies
        no warnings 'prototype';
        no warnings  'redefine';
        *nrg_of_move = sub { 0 };       # call dummy instead of nrg_of_move
    }

    $nrg = RNA::energy_of_struct( $seqstr, pt2str( @pt)) unless defined $nrg;
    my $maxnrg = $nrg + $minh;

    say "genNeighbors: maxnrg=$maxnrg, minh=$minh, noLP=$noLP, shift=$shift" if $verb;
    my @foundPt  = ();              # record finds here
    my @foundNrg = ();              # entry i: energy of struct i in foundPt

    for( my $i=0; $i<@pt; $i++)     # generate all neighbors -- bp indel
    {
        my $j = $pt[$i];
        if( $i < $j)                # paired -> remove bp
        {
            my ($ilg, $olg) = (insGrowLonely(\@pt, $i, $j), outGrowLonely(\@pt, $i, $j));
            if( $noLP && ($ilg || $olg))    # handle lonely pairs (some bp grows lonely)
            {
                    if( isOutLonely(\@pt, $i, $j) && $ilg)    #  bp out-lonely & ins-lonely-growing
                    {
                        my $newnrg = $nrg + nrg_of_move( $seqstr, @pt, $i, $j);
                        remBp( \@pt, $i, $j);
                        $newnrg = $newnrg + nrg_of_move( $seqstr, @pt, $i+1, $j-1);
                        remBp( \@pt, $i+1, $j-1);        # also rem lonely-growing bp
                        addStrIfBetter ($seqstr, \@pt, \@foundPt, \@foundNrg, $maxnrg,
                            { print => $verb ? "     Found delLP: " : "", nrg => $newnrg});
                        addBp( \@pt, $i+1, $j-1);
                        addBp( \@pt, $i, $j);        # Revert changes
                    }
            }
            else    # don't care about lonely pairs (or there is none)
            {
                my $newnrg = $nrg + nrg_of_move( $seqstr, @pt, $i, $j);
                remBp( \@pt, $i, $j);
                addStrIfBetter ($seqstr, \@pt, \@foundPt, \@foundNrg, $maxnrg,
                    { print => $verb ? "     Found (del): " : "", nrg => $newnrg});
                addBp( \@pt, $i, $j);        # Revert changes
            }
        }
        elsif ( $j == -1)     # unpaired -> check all possible pairs
        {
            for my $j (getValBpRight(\@seq, \@pt, $i))  # get valid base pairs
            {
                if( $noLP && isLonely(\@pt, $i, $j) ) # inserted bp lonely=>insert lonely stack
                {
                    if( $j-$i>5 && $pt[$i+1]==-1 && $pt[$j-1]==-1 && isValBp(\@seq,$i+1,$j-1)
                        && isInsLonely( \@pt, $i+1, $j-1))
                    {
                        my $newnrg = $nrg + nrg_of_move( $seqstr, @pt, $i, $j);
                        addBp( \@pt, $i, $j);
                        $newnrg = $newnrg + nrg_of_move( $seqstr, @pt, $i+1, $j-1);
                        addBp( \@pt, $i+1, $j-1);        # also rem lonely-growing bp
                        addStrIfBetter ($seqstr, \@pt, \@foundPt, \@foundNrg, $maxnrg,
                            { print => $verb ? "     Found addLP: " : "", nrg => $newnrg});
                        remBp( \@pt, $i+1, $j-1);        # also rem lonely-growing bp
                        remBp( \@pt, $i, $j);      # Revert changes
                    }    #else: cannot insert stack, skip bp
                }
                else
                {
                    my $newnrg = $nrg + nrg_of_move( $seqstr, @pt, $i, $j);
                    addBp( \@pt, $i, $j);
                    addStrIfBetter ($seqstr, \@pt, \@foundPt, \@foundNrg, $maxnrg,
                        { print => $verb ? "     Found (add): " : "", nrg => $newnrg});
                    remBp( \@pt, $i, $j);      # Revert changes
                }
            }
        }
    }
    for( my $i=0; $shift && $i<@pt; $i++)    # Generate all neighbors -- shift moves
    {
        my $j = $pt[$i];
        next unless  $i < $j;        # paired, i is left end
        my @valbp=();

        if( $noLP ) #Optimise! Only test four cases.
        {
            my $newnrg;
            # Outside: shift i to left or j to right
            if( $pt[$i+1]==$j-1 && !isInsLonely(\@pt, $i+1, $j-1))
            {
                if( $j < $#pt)
                {
                    my $k = $pt[$j+1];    # i -- check only two possible pair (k+1,j)/(j/k+1) [cross]
                    #say "k=$k i=$i k>=0:", $k>=0 , ' k<$#pt:', $k<$#pt, ' k!=i-1:', $k!=$i-1, ' pt[k+1]<0:', $pt[$k+1]<0, ' ValBP(k+1,j):', isValBp(\@seq, $k+1,$j);
                    if( $k>=0 && $k<$#pt && $k!=$i-1 && $pt[$k+1]<0 && isValBp(\@seq, $k+1, $j))
                    {
                        $newnrg = $nrg + nrg_of_move( $seqstr, @pt, $i, $j);
                        remBp( \@pt, $i, $j);
                        #say "Shift i to left / cross: found candidate";
                        my ($lft, $rgt) = $k<$j ? ($k+1, $j) : ($j, $k+1);
                        my $newernrg = $newnrg + nrg_of_move( $seqstr, @pt, $lft, $rgt);
                        addBp( \@pt, $k+1, $j);
                        addStrIfBetter ($seqstr, \@pt, \@foundPt, \@foundNrg, $maxnrg,
                            { print => $verb ? "     Found silLP: " : "", nrg => $newernrg});
                        remBp( \@pt, $k+1, $j);
                        addBp( \@pt, $i, $j);
                    }
                }
                if( $i > 0)
                {
                    my $k = $pt[$i-1];    # j -- check only two possible pair (i,k-1)/(k-1,i) [cross]
                    if( $k>0 && $k!=$j+1 && $pt[$k-1]<0 && isValBp(\@seq, $i, $k-1))
                    {
                        $newnrg = $nrg + nrg_of_move( $seqstr, @pt, $i, $j) unless defined $newnrg;
                        remBp( \@pt, $i, $j);
                        #say "Shift j to right / cross: found candidate";
                        my ($lft, $rgt) = $i<$k ? ($i, $k-1) : ($k-1, $i);
                        my $newernrg = $newnrg + nrg_of_move( $seqstr, @pt, $lft, $rgt);
                        addBp( \@pt, $i, $k-1);
                        addStrIfBetter ($seqstr, \@pt, \@foundPt, \@foundNrg, $maxnrg,
                            { print => $verb ? "     Found sjrLP: " : "", nrg => $newernrg});
                        remBp( \@pt, $i, $k-1);
                        addBp( \@pt, $i, $j);
                    }
                }
            }
            # Inside: shift i to right or j to left
            if( $i>0 && $pt[$i-1]==$j+1 && !isOutLonely(\@pt, $i-1, $j+1))
            {
                my $k = $pt[$j-1];    # i -- check only possible pair (k-1,j)
                if( $k>$i+1 && $pt[$k-1]<0 && isValBp(\@seq, $k-1, $j))
                {
                    $newnrg = $nrg + nrg_of_move( $seqstr, @pt, $i, $j) unless defined $newnrg;
                    remBp( \@pt, $i, $j);
                    #say "Shift i to right: found candidate";
                    my $newernrg = $newnrg + nrg_of_move( $seqstr, @pt, $k-1, $j);
                    addBp( \@pt, $k-1, $j);
                    addStrIfBetter ($seqstr, \@pt, \@foundPt, \@foundNrg, $maxnrg,
                        { print => $verb ? "     Found sirLP: " : "", nrg => $newernrg});
                    remBp( \@pt, $k-1, $j);
                    addBp( \@pt, $i, $j);
                }
                $k = $pt[$i+1];    # j -- check only possible pair (i,k+1)
                if( $k>=0 && $k<$j-1 && $pt[$k+1]<0 && isValBp(\@seq, $i, $k+1))
                {
                    $newnrg = $nrg + nrg_of_move( $seqstr, @pt, $i, $j) unless defined $newnrg;
                    remBp( \@pt, $i, $j);
                    #say "Shift j to left: found candidate";
                    my $newernrg = $newnrg + nrg_of_move( $seqstr, @pt, $i, $k+1);
                    addBp( \@pt, $i, $k+1);
                    addStrIfBetter ($seqstr, \@pt, \@foundPt, \@foundNrg, $maxnrg,
                        { print => $verb ? "     Found sjlLP: " : "", nrg => $newernrg});
                    remBp( \@pt, $i, $k+1);
                    addBp( \@pt, $i, $j);
                }
            }
        }
        else
        {
            my $newnrg = $nrg + nrg_of_move( $seqstr, @pt, $i, $j);
            remBp( \@pt, $i, $j);
            # Shift i
            @valbp = getValBpLeft(\@seq, \@pt, $j);
            @valbp = grep { $_ != $i } @valbp;        # skip (i,j)
            push @valbp, getValBpRight(\@seq, \@pt, $j); # slide i beyond j, now j<i!
            #say "($i,$j): Found ", scalar(@valbp), " valid slide moves for i:\n\t@valbp";
            for my $k (@valbp)     # slide i
            {
                #say "($i,$j): shift $i to $k, lonely=",  isLonely(\@pt, $k, $j),
                #    " insLone=", isInsLonely(\@pt,$k,$j), " outLone=", isOutLonely(\@pt,$k,$j);
                my ($lft, $rgt) = $k<$j ? ($k, $j) : ($j, $k);
                my $newernrg = $newnrg + nrg_of_move( $seqstr, @pt, $lft, $rgt);
                addBp( \@pt, $k, $j);
                addStrIfBetter ($seqstr, \@pt, \@foundPt, \@foundNrg, $maxnrg,
                    { print => $verb ? "     Found (shi): " : "", nrg => $newernrg});
                remBp( \@pt, $k, $j);
            }

            # Shift j
            @valbp = getValBpRight(\@seq, \@pt, $i);
            @valbp = grep { $_ != $j } @valbp;        # skip (i,j)
            push @valbp, getValBpLeft(\@seq, \@pt, $i); # slide j beyond i, now j<i!
            #say "($i,$j): Found ", scalar(@valbp), " valid slide moves for j:\n\t@valbp";
            for my $k (@valbp)     # slide j
            {
                #say "($i,$j): shift $j to $k, lonely=",  isLonely(\@pt, $i, $k),
                #    " insLone=", isInsLonely(\@pt,$k,$i), " outLone=", isOutLonely(\@pt,$k,$i);
                my ($lft, $rgt) = $i<$k ? ($i, $k) : ($k, $i);
                my $newernrg = $newnrg + nrg_of_move( $seqstr, @pt, $i, $k);
                addBp( \@pt, $i, $k);
                addStrIfBetter ($seqstr, \@pt, \@foundPt, \@foundNrg, $maxnrg,
                    { print => $verb ? "     Found (shj): " : "", nrg => $newernrg});
                remBp( \@pt, $i, $k);
            }
            addBp( \@pt, $i, $j);
        }

    }
    {   # Restore function
        no warnings 'prototype';
        no warnings  'redefine';
        *nrg_of_move = $nrg_of_move_sub unless $do_energies;
    }
    return () unless @foundNrg;    # No matching neighbors found, return empty list ()

    parSort(\@foundNrg, \@foundPt) if $sort;
    my @found=();
    for my $i( 0..$#foundNrg)
    {
        push @found, $foundPt[$i], $foundNrg[$i];
    }
    return @found;
}


# Wrapper for genNeighbors function. genNeighbors returns a list of form
# (struct1, energy1, struct2, energy2, ...) which is performant, but
# impractical as one cannot easily iterate over it non-destructively (i.e.
# without while splice x, 0, 2). This sub wraps the list up into pairs
# ([struct1, energy1], [struct2, energy2], ...)
# Arguments: passed on to genNeighbors.
sub genNeighbors_pairs {
    return pairs genNeighbors @_;
}


# Return 1 iff passed bp is inside-lonely, i.e. has no inner neighbors.
sub is_ins_lonely ( +$$ ) { isInsLonely(@_) }       # snake case wrapper
sub isInsLonely ( +$$ )
{
    my ($ptref, $i, $j) = @_;
    # $j = $ptref->[$i] unless defined($j);
    assert defined $j;
    #($i,$j)=($j,$i) if $j<$i;
    assert $i<$j;
    return $ptref->[$i+1] != $j-1;
}


# Return 1 iff passed bp is outside-lonely, i.e. has no surrounding neighbor.
sub is_out_lonely ( +$$ ) { isOutLonely(@_) }       # snake case wrapper
sub isOutLonely ( +$$ )
{
    my ($ptref, $i, $j) = @_;
    #$j = $ptref->[$i] unless defined($j);
    assert defined $j;
    #($i,$j)=($j,$i) if $j<$i;
    assert $i<$j;
    return $i==0 || $j==$#$ptref || $ptref->[$i-1] != $j+1;
}


# Returns true iff position i is paired to another position.
sub is_paired ($$) {
    my ($pt, $i) = @_;
    my $is_paired = $pt->[$i] >= 0;
    return $is_paired;
}

# Check whether base pair (i, j) in a given pair table is lonely.
sub is_lonely ( +$$ ) { isLonely(@_) }      # snake case wrapper
sub isLonely ( +$$ )
{
    my ($ptref, $i, $j) = @_;
    return isOutLonely( $ptref, $i, $j) && isInsLonely( $ptref, $i, $j);
}


# Check whether a given RNA secondary structure is canonical, i.e. if it does
# not contain any lonely base pairs (cf. isLonely). The structure may be given
# as a dot-bracket string or a 0-based pair table.
# Returns true if structure is canonical, and false if it is not.
sub is_lonely_struct ($) { isCanonical(@_) }
sub is_canonical         { isCanonical(@_) }            # snake case wrapper
sub isCanonical {
    my ($struct_or_pt) = @_;
    my $pt = reftype $struct_or_pt ? $struct_or_pt : str2pt $struct_or_pt;

    my $is_canonical
        = none {my $j = $pt->[$_]; $j>$_ and isLonely $pt, $_, $j} 0..$#$pt;
    return $is_canonical;
}

# TODO deprecated. There is no need for a sequence here ... use
# is_lonely_struct() instead.
sub isLonelyStruct( $$ )
{
    my ($seq, $str) = @_;
    my @pt = str2pt( $str);
    for( my $i=0; $i<@pt; $i++)
    {
        my $j = $pt[$i];
        return 1 if $j>=0 && $i<$j && isLonely( \@pt, $i, $j);
    }
    return 0;
}


# Determine whether there is a bp enclosing (i,j) which will become
# a lonely pair ("grow lonely") when deleting (i,j).
# Returns 0 if no such pair exists; 1 else.
sub outGrowLonely( +$$ )
{
    my ($ptref, $i, $j) = @_;
    return 1 if $i!=0 && $ptref->[$i-1]==$j+1 && isOutLonely( $ptref,$i-1,$j+1);
    return 0;
}


# Determine whether there is a bp enclosed by (i,j) which will become
# a lonely pair ("grow lonely") when deleting (i,j).
# Returns 0 if no such pair exists; 1 else.
sub insGrowLonely( +$$ )
{
    my ($ptref, $i, $j) = @_;
    return 1 if $ptref->[$i+1]==$j-1 && isInsLonely( $ptref,$i+1,$j-1);
    return 0;
}


# Determine whether there is a bp enclosing or enclosed by (i,j) which
# will become a lonely pair ("grow lonely") when deleting (i,j).
# Returns 0 if no such pair exists; 1 else.
sub growLonely( +$$ )
{
    my ($ptref, $i, $j) = @_;
    return 1 if insGrowLonely($ptref,$i,$j) or outGrowLonely($ptref,$i,$j);
    return 0;
}


# Add bp (i,j) to passed passed pair table.
sub addBp( +$$)
{
    my ($ptref, $i, $j) = @_;
    assert $ptref->[$i]==-1 || $ptref->[$j]==-1; # at least one unpaired, e.g. sliding moves
    $ptref->[$i] = $j;
    $ptref->[$j] = $i;
}


# Add bp (i,j) to passed passed pair table.
sub remBp( +$$)
{
    my ($ptref, $i, $j) = @_;
    #$j=$ptref->[$i] unless defined($j);
    assert defined $j;
    assert $ptref->[$i]>=0 || $ptref->[$j]>=0; # at least one paired, e.g. sliding moves
    $ptref->[$i] = -1;
    $ptref->[$j] = -1;
}


# Add structure and its energy to passed lists if better than passed max.
# If optional argument 'print' is specified and non-empty, print the
# message from 'print', structure and energy if structure better than max.
# strnrg: optional, energy of passed structure. Computed if not passed
# opts: hash reference for optional arguments
# Return 1 if structure was added / is better and 0 else.
sub addStrIfBetter
{
    my ($seqstr, $ptref, $foundPtRef, $foundNrgRef, $maxnrg, $opts) = @_;
    my ($newnrg, $print) = @{$opts}{qw(nrg print)};    # hashref slice => shitty syntax
    $newnrg = RNA::energy_of_struct( $seqstr, pt2str( @$ptref) ) unless defined $newnrg;
    #say "   Next candidate      : ", printStrNrg(\@pt, $newnrg);
    if( $newnrg <= $maxnrg)
    {
        say "$print", printStrNrg($ptref, $newnrg) if $print;
        my @newpt = @$ptref;    # Copy pairtable
        push @$foundPtRef, \@newpt;
        push @$foundNrgRef, $newnrg;
        return 1;
    }
    return 0;
}


# Return 1 if (i,j) is valid base pair, i.e. distance of i and j
# is >4 (3 interior bases) and the bases at i and j can pair; 0 else.
sub isValBp
{
    my ($seqref, $i,$j) = @_;
    #say "abs($j-$i)>3:", abs($j-$i)>3, '"$seqref->[$i]$seqref->[$j]"~~@valbp:', "$seqref->[$i]$seqref->[$j]"~~@valbp;
    return 1 if abs($j-$i)>3
                and any {"$seqref->[$i]$seqref->[$j]" eq $_} @VALID_BP;
    return 0;
}


# Given the sequence, pair table and starting index i, generate a list
# of all j such that (i,j) is a valid base pair with i<j
sub getValBpRight
{
    my ($seqref, $ptref, $i) = @_;
    assert $ptref->[$i]<0;
    assert scalar(@VALID_BP);
    my @bp = ();                       # store found bp here
    for( my $j=$i+1; $j<@$ptref; $j++) # start next to $i to count all brackets
    {
        my $k = $ptref->[$j];
        if( $j < $k)         # j is an opening bracket
        {
            $j = $k;
            next;
        }
        elsif( $k>=0)       # j is a closing bracket, end of current loop!
        {
            last;
        }
        push @bp, $j if isValBp($seqref,$i,$j);
    }
    return @bp;
}


# Given the sequence, pair table and starting index j, generate a list
# of all i such that (i,j) is a valid base pair with i<j.
sub getValBpLeft
{
    my ($seqref, $ptref, $j) = @_;
    assert $ptref->[$j]<0;
    my @bp = ();                       # store found bp here
    my $nbrack = 0;
    for( my $i=$j-1; $i>=0; $i--) # start next to $i to count all brackets
    {
        if( $i < $ptref->[$i])         # i is an opening bracket
        {
            $nbrack++;
            last if $nbrack>0;         # reached start of enclosing loop, end
            next;
        }
        elsif( $ptref->[$i]!=-1)       # i is a closing bracket
        {
            $nbrack--;
            next;
        }
        # minimal loop length is 3; crossing pair; invalid bp
        next if $nbrack || !isValBp($seqref,$i,$j);    # Crossing pair?
        push @bp, $i;
    }
    return @bp;
}

# Makes an iterator for a given a pair table, which returns, at each call, the
# next base pair (i,j) contained in the structure encoded by the pair table.
# The base pairs are sorted by the first index.
# Returns next bp (i, j) or false after the last base pair has been returned.
sub bp_iter_factory {
    my ($pt_ref) = @_;
    my $i = -1;
    my $seq_len = seq_len $pt_ref;

    return sub {
        while ($i < $seq_len-1) {           # run till last index after inc
            $i += 1;                        # prepare for next call

            # If pt[i] >= 0, i is paired to some j. Also exclude pair (j, i)
            # where j > i.
            my $j = $pt_ref->[$i];
            next unless $i < $j;

            return ($i, $j);
        }
        return;                             # no more base pairs
    };
}

1;     # Required to return 1 from a Perl package...
# EOF




__END__

=head1 NAME

RNAhelper - Utility and misc functions for RNA structure bioinformatics.

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module is a collection of misc functions for repeating tasks in RNA
structure bioinformatics, e. g. mfe folding, dot--bracket string parsing,
gradient walks on structures, Boltzmann weight computation etc. pp.

Here is an example:

    use v5.12;
    use RNAhelper qw(rna_mfe_struct);

    say rna_mfe_struct 'GGGAUGCCC';


=head1 EXPORTS

Functions need to be imported explicitly. This is a list of all exported
functions. A more detailed description is found in the next section.

=over

=item printStrNrg

=item printStrNrgLst

=item pt2str

=item pt2str_c

=item str2pt

=item path2str

=back

=over

=item rand_seq

=item rand_seed

=back

=over

=item grad_walk

=item floodCycle

=back

=over

=item bestDirectPath

=item bestDirectSaddle

=item directPaths

=back

=over

=item genNeighbors

=item genNeighbors_pairs

=item isLonely, is_lonely

=item isLonelyStruct, is_lonely_struct

=item isCanonical, is_canonical

=item isInsLonely, is_ins_lonely

=item isOutLonely, is_out_lonely

=back

=over

=item is_valid_seq

=item normalize_seq

=item seq_len

=back

=over

=item bpValid, bp_valid

=item is_paired

=item struct_cmp

=item struct_lt

=item struct_eq

=item bp_dist

=item bp_iter_factory

=item nrg_of_move

=back

=over

=item rna_mfe

=item rna_mfe_struct

=item rna_eval_structure

=item make_fold_compound

=item make_model_details

=back

=over

=item boltz2nrg

=item nrg2boltz

=item ensemble_energy

=item subopt_energy

=item partition_energy

=item partition_energy_structiter

=item partition_energy_nrgiter

=item get_gas_const

=item celsius2kelvin

=item get_basin_partition_funcs

=back

=over

=item base_pair_prob

=back

=over

=item get_opt

=item parUniq

=item parSort

=back

=over

=item test

=back


=head1 SUBROUTINES

Detailed documentation of exported subroutines.

=head2 printStrNrg($ptref | $dbstring, $nrg)

Convert pairtable and energy value to a formatted string. Pass structure and
energy of structure. Excepts either a reference to a pairtable or a
dot-bracket string as first arg.

=head2 printStrNrgLst($pt_energy_list, $prefix = "")

Print a list containing (pair table ref, energy) pairs,
each entry prefixed with an optionally passed prefix.

=head2 pt2str, pt2str_c, struct_string($pt_ref | @pair_table)

Converts a pair table to a dot-bracket string.

Arguments:

=over

=item pair_table: list or array ref to pairing table.

=back

Returns dot-bracket string of the structure.

=head2 str2pt($dot_bracket_string)

Generate pair table from dot-bracket string.
Indexing starts at 0, -1 means unpaired (unlike ViennaRNA!). Returns an
arrayref in scalar context, or a list of base pairs in list context.

=head2 path2str($path)

For the passed path, return a list of dot-bracket strings of visited structures
from start to target.

=head2 rand_seq($length, $gc_content = undef)

Generate a random RNA sequence consisting of G, C, A and U. Optionally, the
exact GC content of the sequence can be specified. If not specified, it is
chosen randomly.

Arguments:

=over

=item length: ~ of the sequence to be generated.

=back

Optional arguments:

=over

=item gc_content: Expected GC content of the generated sequence. 'Expected'
    means that this value is used as a probability for sampling either a G
    or a C, so the actual GC content of the generated sequences may vary

=back

Returns a random sequence with the requested GC content.

=head2 rand_seed($seed = undef)

Seed the random number generator (Mersenne Twister) with the passed value.
This is ONLY required to reproduce results, not before normal use.
When passing undef, the generator is auto-seeded with a random integer. The
used seed is returned. This is FAR LESS random than the random generator's
internal auto-seeding, so really don't use this unless you need to.

Arugments:

=over

=item seed: The value to be passed to Math::Random::MT::Auto->srand. Can be a
      single integer, an array ref containing integers etc. Passing undef
      will auto-seed the generator with a single integer.

=back

Returns the used seed.

=head2 grad_walk

Perform a gradient walk from current structure down to local minimum of the
current basin. Correctly handle degenerate landscapes, cf. grad_step

Arguments:


=over

=item sequence: sequence string

=item start_struct:  reference to pairtable of start structure

=item start_energy:  energy of start structure (optional, is computed if not
               defined)

=item opt: option hash ref:

=over

=item walk_structs: array ref used to store structs (pt's) and energies of
    structures encountered on walk between first and last struct if
    they are required

=item verb: flag, be verbose? [0]

=item noLP: flag, use moveset avoiding lonely pairs? [0]

=item shift: flag, use moveset using shift moves? [1]

=back

=back

Returns the last structure (pairtable) of the walk and its energy.

=head2 floodCycle($seqstr, $startptref, $minh, $qedref, $minNmin, $maxFlood, $ofh, $noLP)

Cycle (basin) exploration, algorithm description in summary_riboswitches!

Flood cycle up to minh above starting structure's minimum and report all
surrounding local minima reachable from transient structures found.

The following optional params are mainly used by waterfall:

=over

=item qedref:   optional, pass ref to global hash of visited structures to re-use
                this information across multiple calls ("ref to queued hash@).

=item minNmin:  Optional, used with qedref. Only update qedref hash if at least
                <minNmin> minima were found (allows re-visiting structures after
                another minh has been passed.

=item maxFlood: Optional, used with qedref and minNmin. If minh >= maxFlood,
                update qedref in any case, even if not enough minima were found.

=item ofh:      Output file handle, write structures and their energies found in
                this cycle to a file.

=item noLP      Use moveset avoiding lonely pairs by removing/adding lonely
                stacks of size 2 if a single pair indel would be lonely. [1]

=back

=head2 bestDirectPath

Get the best (i.e. energetically lowest) direct path connecting
the start structure with the target structure. The search is heuristic
and keeps at most n structures at a time. The search space can be
restricted to canonical paths (noLP, no lonely pairs). Shift moves
are also supported.

Arguments:

=over

=item seq: Sequence string.

=item start: Start structure as dot-bracket string.

=item target: Target structure as dot-bracket string.

=item n: Keept at most n structures during search.

=item noLP: Do not allow lonely pairs, i.e. search canonical paths [0]

=item shift: Allow shift moves.  [0]

=back

Returns a list of structure strings, including their energies, denoting the
found optimal path.

=head2 bestDirectSaddle

Get the best (i.e. energetically lowest) direct path connecting
the start structure with the target structure and return its saddle
energy. The search is heuristic and keeps at most n structures at a
time. The search space can be restricted to canonical paths
(noLP, no lonely pairs). Shift moves are also supported.

Arguments:

=over

=item seq: Sequence string.

=item start: Start structure as dot-bracket string.

=item target: Target structure as dot-bracket string.

=item n: Keept at most n structures during search.

=item noLP: Do not allow lonely pairs, i.e. search canonical paths.

=item shift: Allow shift moves.

=back

Returns the saddle (i.e. highest) energy of the best (i.e. energetically
lowest) path found, rounded to two digits.

=head2 directPaths

Compute direct paths with low saddles from structure pta to ptb (pass pairtables).
Keep at most n entries during search.

Arguments:

=over

=item seqstr: Sequence as string.

=item pta: Pair table of start structure.

=item ptb: Pair table of target structure.

=item n: Keep at most n structures during search.

=item opt: Optional hash for options, see below.

=back

Options:

=over

=item shift: Allow shift moves. [1]

=item noLP: Only search canonical paths. [1]

=item verb: Be verbose.

=item nreturn: Return only the best <nreturn> paths. [Inf]

=item noDupes: No duplicates. [0]

=item maxNrg: Maximum allowed energy in the path. [Inf]

=back

Returns a list of the paths found.

=head2 genNeighbors

Generate all neighbors and filter w.r.t. their energy difference to the
current structure.

Arguments:

=over

=item seq: Sequence string.

=item pair table ref: Pair table of structure.

=back

Optional arguments:

=over

=item minh:  Maximal energy difference of neighbors. For the energy E of the
             current structure, only return structures with energy at most
             E + minh. [Infinity]

=item noLP:  Skip neighbors containing lonely pairs? [1]

=item shift: Perform shift moves? [1]

=item nrg:   Energy of current structure. Re-computed unless provided.

=item sort:  Sort the returned structure list by energy. [1]

=item do_energies: Compute energy of neighbors? If not, all energy values are
                   set to 0 and options minh, nrg and sort are disabled. [1]

=back

Returns (sorted?) list of pairs (\@pt,$nrg) (ref to pair table, energy)
resembling the neighbor structures fulfilling the energy criterion.

=head2 genNeighbors_pairs

Wrapper for genNeighbors function. genNeighbors returns a list of form
(struct1, energy1, struct2, energy2, ...) which is performant, but
impractical as one cannot easily iterate over it non-destructively (i.e.
without while splice x, 0, 2). This sub wraps the list up into pairs
([struct1, energy1], [struct2, energy2], ...)

Arguments: passed on to genNeighbors.

=head2 isLonely, is_lonely($ptref, $i, $j)

Check whether base pair (i, j) in a given pair table is lonely.

=head2 isLonelyStruct($seq, $str)

Deprecated. Use C<is_lonely_struct()> instead.

=head2 isCanonical, is_canonical, is_lonely_struct($struct_or_pt)

Check whether a given RNA secondary structure is canonical, i.e. if it does
not contain any lonely base pairs (cf. isLonely). The structure may be given
as a dot-bracket string or a 0-based pair table.

Returns true if structure is canonical, and false if it is not.

=head2 isInsLonely, is_ins_lonely($ptref, $i, $j)

Return 1 iff passed bp is inside-lonely, i.e. has no inner neighbors.

=head2 isOutLonely, is_out_lonely($ptref, $i, $j)

Return 1 iff passed bp is outside-lonely, i.e. has no surrounding neighbor.

=head2 is_valid_seq($sequence)

Return true iff a valid sequence is given, that is one consisting solely of
A, U, T, G, and C, either in upper or lower case. The empty sequence is not
valid.

=head2 normalize_seq($sequence)

Normalize sequence, i.e. convert it to upper case, and delete any
white-space characters.

Arguments:

=over

=item sequence: Sequence to normalize.

=back

Returns the normalized sequence.

=head2 seq_len($seq_or_ptref)

Given a sequence string, structure string, or a pair table ref, return the
sequence length.

=head2 bpValid, bp_valid

Check base pair validity. Valid base pairs are AU, UA, GC, CG, GU, and UG.
Pass sequence string and pair table reference.

WARNING: This is slow and should be improved, cf. code.

Returns 1 iff structure is valid.

=head2 is_paired($pair_table, $i)

Returns true iff position i is paired to another position.

=head2 struct_cmp($pt1, $pt2)

Comparison function for two structure pair tables. Let struct1, struct2 be
the two input structures. We define struct1 < struct2 iff
sequence_length(struct1) < sequence_length(struct2) or
sequence_length(struct1) == sequence_length(struct2) and
    struct1 cmp struct2 < 0 (i.e. struct1 is lexicographically smaller)
and use the total order induced by this relation.

Arguments:

=over

=item pt1: Pair table of first structure.

=item pt2: Pair table of second structure.

=back

Returns a comparison value, i.e. a negative value if pt1 < pt2, a positive
value if pt1 > pt2, and zero if pt1 == pt2.

=head2 struct_lt($pt1, $pt2)

Check whether a structure is smaller than another one. For details, cf.
struct_cmp.

Arguments:

=over

=item pt1: (Ref to) pair table of first structure.

=item pt2: (Ref to) pair table of second structure.

=back

Returns true value if the first structure is smaller than the second one,
and a false value otherwise.

=head2 struct_eq

Check whether two structures are equal. For details, cf. struct_cmp.

Arguments:

=over

=item pt1: (Ref to) pair table of first structure.

=item pt2: (Ref to) pair table of second structure.

=back

Returns true value if the structures are equal, and a false value otherwise.

=head2 bp_dist($struct_a, $struct_b)

Compute the base pair distance of two structures given either as pair tables
or dot-bracket strings.

=head2 bp_iter_factory($pt_ref)

Makes an iterator for a given a pair table, which returns, at each call, the
next base pair (i,j) contained in the structure encoded by the pair table.
The base pairs are sorted by the first index.
Returns next bp (i, j) or false after the last base pair has been returned.

=head2 nrg_of_move($seqstr, $ptref, $i, $j)

Wrapper for RNA::energy_of_move from ViennaRNA to efficiently compute the
energy difference of a single bp move.

Arguments:

=over

=item seqstr: sequence string

=item ptref:  reference to pairtable of struct

=item i, j:   base pair that WILL BE BUT IS NOT YET added/removed to/from
              pairtable

=back

Returns energy difference that would arise when adding / removing this bp.

Notes:

=over

=item ViennaRNA's pairtables are 1-based, in this package they are 0-based.

=item ViennaRNA tries to ADD (i,j) if i,j>0, and to REMOVE (i,j) if i,j<0.

=item ViennaRNA requires i<j or the package could die

=back

=head2 rna_mfe($seq, $opt_ref)

Compute the minimum free energy (mfe) of a given sequence, considering
optional arguments like temperature and dangling end model.

Arguments:

=over

=item seq: RNA sequence to compute mfe for.

=item opt_ref: Hash ref containing options, see make_model_details().

=back

In scalar context, returns the mfe of the sequence. In list context, both
the mfe structure and the MFE itself are returned.

=head2 rna_mfe_struct($seq, $opt_ref)

Compute the structure of minimum free energy for a given sequence,
considering optional arguments like temperature and dangling end model.

Arguments:

=over

=item seq: RNA sequence to compute mfe structure for.

=item opt_ref: Hash ref containing options, see make_model_details().

=back

Returns the mfe structure as dot-bracket string.

=head2 rna_eval_structure($seq, $struct, $opt_ref)

Compute the free energy of an RNA secondary structure for a given sequence,
considering optional arguments like temperature and dangling end model.

Arguments:

=over

=item seq: RNA sequence to compute energy for

=item struct: secondary structure to evaluate

=item opt_ref: hash ref containing options, see make_model_details().

=back

Returns the energy of the given structure.

=head2 make_fold_compound($seq, $opt_ref)

Create a ViennaRNA fold compound object from a sequence. If model options
are given as a hash ref, create a model detail object and integrate settings
into the fold compound.

Arguments:

=over

=item seq: RNA sequence to create compound for.

=back

Optional args, given as hash ref: cf. make_model_details().

Returns ViennaRNA fold compound of the given sequence and settings.

=head2 make_model_details($opt_ref)

Create a model details object which can be passed to ViennaRNA's fold
compound constructor to set model details for the energy computation.

Arguments: hash ref containing:

=over

=item temperature: the RNA folding temperature in deg C [37]

=item dangles: dangling end model, either 0 (= no dangles), 1 , 2
      (full/overdangle), or 3 (+ helix stacking) [2]

=item noLP: prohibit lonely base pairs (1 - only canonical structs, 0 - all) [0]
            DOES NOT WORK AS EXPECTED for partition function calculations!
            E.g. subopts will still contain some structures with lonely pairs,
            and partition function will also not be exact.

=item pf_smooth: smooth energy parameters for partition function calculations [1]
                 Important to have differentiable energies when varying
                 temperature. Disable to get exact partition functions when
                 comparing with individual structures' Boltzmann weight (e.g.
                 for coverage analyses, probabilities of structures etc.)

=item compute_bpp: Compute base pair probabilities when doing partition function
                   folding. [1] Disable improve performance.

=back

Returns the model details object.

=head2 boltz2nrg($boltzmann_weight, $temperature)

Compute the energy value E associated with the given Boltzmann weight or
partition function Z, i.e.

      E = - R * T * ln Z

where T is the temperature and R is the universal gas constant.

Arguments:

=over

=item boltzmann_weight: Boltzmann weight or partition function for which to
          compute the associated energy value

=back

Optional arguments:

=over

=item temperature: temperature at which to compute the Boltzman weight, given in
          degree Celsius

=back

Returns the energy value.

=head2 nrg2boltz($energy, $temperature)

Compute the Boltzmann weight Z of a given energy value E, i.e.

      Z = exp( - E / (R * T) ),

where T is the temperature and R is the universal gas constant.

Arguments:

=over

=item energy: energy of which to compute the Boltzmann weight, given in kcal/mol

=back

Optional arguments:

=over

=item temperature: temperature at which to compute the Boltzman weight, given in
          degree Celsius

=back

Returns the Boltzmann weight / factor.

=head2 ensemble_energy($seq, $opt_ref)

Compute ensemble / partition free energy. It does the same as
partition_energy($seq), but accepts RNA model detail options to adjust
temperature or dangling model.

NOTE: Use parameter pf_smooth => 0 when using the ensemble energy together
with energies of single structures (e.g. when computing probabilities).
Smoothing would lead to errors!

Arguments:

=over

=item seq:         sequence for which the ensemble energy should be computed

=back

Optional args passed in hash ref: cf. make_model_details().
The optional arguments default to their ViennaRNA default values.

Returns the (full) ensemble free energy.

=head2 subopt_energy

Compute the partial ensemble energy of an enumerated energy band above the
mfe structure.

Arguments:

=over

=item seq:        Sequence for which the ensemble energy should be computed.

=item e_thresh:   Relative enumeration threshold (w.r.t. mfe) that is supposed
                  to be enumerated.

=back

Optional args passed in hash ref: cf. make_model_details().
The optional arguments default to their ViennaRNA default values.

NOTE: The subopt noLP (no lonely pair) option does NOT WORK THE SAME as the
implementations in this module. Specifically, some structures may still
contain lonely pairs.

Returns the (partial) ensemble free energy of the enumerated energy band.

=head2 partition_energy($seq, $opt_ref)

Wrapper for the partition_energy_iter function that can be called using a
a sequence and an array of structures.

Arguments:

=over

=item seq: Sequence for which to compute the structure energies.

=back

Optional arguments:

=over

=item opt_ref: Model detail options, cf. make_model_details().

=item structures: List of secondary structure dot--bracket strings for which to
          compute the free energy. If empty, computes energy of entire ensemble.

=back

Returns the total energy covered by the passed structures. If none were
passed, return the full ensemble energy.

=head2 partition_energy_structiter($seq, $struct_iter, $scale_energy, $opt_ref)

Given a sequence and a struct iterator, computes the partition energy of the
set of structures returned by the iter. Wrapper for
C<partition_energy_structiter()> that accepts a structure iter instead of an
energy iter.

Arguments:

=over

=item seq: Sequence for which to compute the structure energies.
      struct_iter: Iterator that returns the next structure for which to add up
          the energies.

=back

Optional arguments:

=over

=item scale_energy: Energy offset which is substracted from all structures'
          energies before their Boltzmann weight is calculated. This avoids
          numerical issues. Default: mfe of sequence.

=item opt_ref: hash ref containing additional folding energy options:
          temperature: folding temperature for Boltzmann weight conversions

=back

Returns the scaled partition energy of the given structures.

=head2 partition_energy_nrgiter($seq, $energy_iter, $scale_energy, $opt_ref)

Computes the "partition energy" for a list of free energies. This is the
energy related to the partition function by the Boltzmann transformation.
More precisely, for a set of structures X and its partition function Z given
by Z = sum Boltz( E(x) ) over x in X, where E(x) is the energy of x, its
partition energy E_p and is computed such that Z = Boltz( E_p ).

Arguments:

=over

=item seq: Sequence for the given list of structures.
      energy_iter: iterator that returns, with each call, the next secondary
          structure's energy.

=back

Optional arguments:

=over

=item scale_energy: Energy offset which is substracted from all structures'
          energies before their Boltzmann weight is calculated. This avoids
          numerical issues. Default: mfe of sequence

=item opt_ref: hash ref containing additional folding energy options:
          temperature: folding temperature for Boltzmann weight conversions

=back

Returns the partition energy of the given list of free energies.

=head2 get_gas_const

Returns the universal gas constant in kcal/mol as defined in ViennaRNA.

=head2 celsius2kelvin($deg_celsius)

Convert a temperature value from degree Celsius to Kelvin.

=head2 get_basin_partition_funcs($sbmap_filehandle, $opt)

Compute the partition function of each basin, normalized with the mfe.

Arguments:

=over

=item sbmap_filehandle: Handle to a file that lists on each line a  structure,
        its energy and its basin number; separated by whitespace. The structure
        is not used and can be replaced by a non-whitespace dummy.

=item temperature: Folding temperature in deg Celsius. [37]

=item mfe: Optional, minimum free energy of sequence as a scalar. If undefined,
        it is determined by reading the first line from the sbmap file (assumes
        ordering by energy!). If a reference is provided, the mfe is determined
        from sbmap_filehandle and the referenced scalar is updated with the
        result. Otherwise, the given mfe is used.

=back

Returns a list of scaled basin partition functions, where the index
corresponds to the basin number.

=head2 base_pair_prob

Given an RNA sequence of length n, compute an n x n probability matrix bpp
storing the probabilities that nucleotide i is paired with nucleotide j in
entry bpp[i,j]. This matrix is symmetric.

Instead of a sequence, a ViennaRNA fold compound may be passed to allow
adding constraints, etc. to the sequence.

Arguments:

=over

=item seq_or_fc: Nucleotide sequence or ViennaRNA fold compound containing the
    sequence.

=back

Returns reference to 2-dim array containing probabilities that index i pairs
with index j.

=head2 get_opt

Extract a value from an option hash. If the option does not exist in the
passed hash, a default value is returned if it was passed, otherwise undef.

Arguments:

=over

=item opt_ref:    Reference to a hash containing options.

=item opt_name:   String of the options name.

=item default:    Optional, default value if opt_name is not defined.

=back

Returns the value of the option if it exists in the opt_ref hash, otherwise
a default value if it was passed, otherwise undef.

=head2 parUniq($list_1_ref, ...)

Parallely make passed (SORTED!) lists of equal length unique w.r.t. first list.
More precisely, in the first list, every item is compared to the previous
one, and if both are equal, it is discarded. If position i is discarded in
the first list, it is also discarded in all other lists.

The passed lists are modified in place.

=head2 parSort($list_1_ref, ...)

Parallely sort lists of equal length w.r.t. first list.
More precisely, the first list is sorted numerically, and to the other
lists, the same ordering is applied.

The lists are sorted in place.

=head1 AUTHOR

Felix Kuehnl, C<< <felix at bioinf.uni-leipzig.de> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-rnahelper at rt.cpan.org>,
or through the web interface at
L<https://rt.cpan.org/NoAuth/ReportBug.html?Queue=RNAhelper>.  I will be
notified, and then you'll automatically be notified of progress on your bug as
I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc RNAhelper


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<https://rt.cpan.org/NoAuth/Bugs.html?Dist=RNAhelper>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/RNAhelper>

=item * CPAN Ratings

L<https://cpanratings.perl.org/d/RNAhelper>

=item * Search CPAN

L<https://metacpan.org/release/RNAhelper>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2015--2022 Felix Kuehnl.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see L<http://www.gnu.org/licenses/>.


=cut

1; # End of RNAhelper
