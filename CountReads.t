#!/usr/bin/env perl

use strict;
use warnings;

use Test::More;
use File::Temp;
use File::chdir;
use File::Basename;
use Cwd 'abs_path';
require_ok('CountReads');


my $header_line = "reference,READ,G241T,C290T";
my $positions = [ 
    { variant => 'T', id => 'G241T' },
    { variant => 'T', id => 'C290T' },
    ];

my $parsed_variants = CountReads->parse_header($header_line);
is_deeply($parsed_variants, $positions, "Header parsed correctly");

ok(CountReads->valid(qw( G C )), "Valid read identified");
ok(!CountReads->valid(undef, undef), "Non-spanning read invalid");
ok(CountReads->valid(qw( G N )), "Read with N is valid");

is(CountReads->variant_for_base($positions->[0], 'T'), $positions->[0], "Variant identified correctly");
isnt(CountReads->variant_for_base($positions->[0], 'G'), $positions->[0], "Non-variant identified correctly");

my @variants = CountReads->variants_in_read($positions, 'G', 'T');
is_deeply(\@variants, [$positions->[1]], "Variants returned correctly");

is(CountReads->create_readname_filename("G16T,G30A,C253A"), "G16T_G30A_C253A_readnames.txt", "Filenames constructed correctly");

is(CountReads->correct_readname("read1"), "read1", "Readname correction doesn't screw up test data");
is(CountReads->correct_readname("read1/1333/ccs"), "read1/1333", "Readname correction works");

#integration test
my $src_dir = abs_path(dirname(__FILE__));
my $test_data = "$src_dir/test_data/";
my $input = "$test_data/test_read_variants.csv";
my $expected_output = "$test_data/expected_output.txt";

subtest "CLI Integration Test" => sub {
    my $tempdir = File::Temp->newdir();
    local $CWD = $tempdir;
    system("perl $src_dir/CountReads.pm $input > test_output.csv 2>/dev/null") == 0
        or die "Error running integration test\n";

    ok(!`diff $expected_output test_output.csv`, "Expected output produced");
    for my $file (glob "*_readnames.txt") {
        ok(!`diff $test_data/expected_$file $file`, "Expected readnames match in $file");
    }
};

subtest "Package Integration Test" => sub {
    my $tempdir = File::Temp->newdir();
    local $CWD = $tempdir;

    my $ifh = IO::File->new($input, "r") or die "Unable to open filehandle to read $input\n";
    my $ofh = IO::File->new("test_output2.csv", "w") or die "Unable to open filehandle to write test_output2.csv\n";
    CountReads->run($ifh, $ofh) == 0
        or die "Error running integration test with filenames\n";
    ok(!`diff $expected_output test_output2.csv`, "Expected output produced");
    for my $file (glob "*_readnames.txt") {
        ok(!`diff $test_data/expected_$file $file`, "Expected readnames match in $file");
    }
};

done_testing();
