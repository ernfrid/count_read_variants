#!/usr/bin/env perl

use strict;
use warnings;

use Test::More;

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
ok(!CountReads->valid(qw( G N )), "Read with N is invalid");

is(CountReads->variant_for_base($positions->[0], 'T'), $positions->[0], "Variant identified correctly");
isnt(CountReads->variant_for_base($positions->[0], 'G'), $positions->[0], "Non-variant identified correctly");

my @variants = CountReads->variants_in_read($positions, 'G', 'T');
is_deeply(\@variants, [$positions->[1]], "Variants returned correctly");

done_testing();
