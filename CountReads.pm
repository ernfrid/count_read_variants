#!/usr/bin/env perl
package CountReads;

__PACKAGE__->run() unless caller();

use strict;
use warnings;
use List::MoreUtils qw(pairwise);

sub run() {
    my $class = shift;
    my $variants;
    my %variant_counts;
    my $discarded = 0;
    my $spanning_counts = 0;
    my @position_counts;

    while(<>) {
        $_ =~ s/\R$//g;
        unless($variants) {
            $variants = $class->parse_header($_);
            next;
        }
        my ($ref_name, $read_name, @bases) = split ",", $_, -1;
        if($class->valid(@bases)) {
            my @variant_in_read = $class->variants_in_read($variants, @bases);
            if(@variant_in_read) {
                my $name = join(",", map {$_->{id}} @variant_in_read);
                $variant_counts{$name} += 1;
            }
            $spanning_counts += 1;
        }
        else {
            $discarded += 1;
        }
        for my $index (0..$#bases) {
            $position_counts[$index]->{$bases[$index]} += 1;
        }
    }
    for(my $index = 0; $index < @$variants; ++$index) {
        print STDERR $variants->[$index]->{id}, ":\n";
        for my $base (keys %{$position_counts[$index]}) {
            print STDERR "\t$base\t$position_counts[$index]->{$base}\n";
        }
    }
    print STDERR "Discarded $discarded non-spanning reads\n";
    print STDERR "Evaluated $spanning_counts spanning reads\n";
    print "Variant(s)\tVariant Supporting Reads\tTotal Spanning Reads\tVAF\n";
    for my $variant (sort keys %variant_counts) {
        printf "%s\t%d\t%d\t%f%%\n", $variant, $variant_counts{$variant}, $spanning_counts, $variant_counts{$variant} / $spanning_counts * 100;
    }
}

sub parse_header {
    my ($class, $header_line) = @_;
    chomp $header_line;

    my @parsed_variants;

    my ($ref_name, $read_name_header, @variants) = split ",", $header_line;
    for my $variant (@variants) {
        my ($ref, $pos, $var) = $variant =~ /(\D+)(\d+)(\D+)/;
        unless(defined $ref and defined $var) {
            die "Error parsing variant $variant from header\n";
        }
        push @parsed_variants, {id => $variant, variant => $var};
    }
    return \@parsed_variants;
}

sub valid {
    my ($class, @bases) = @_;
    for my $base (@bases) {
        if( !defined $base
            || $base eq ''
            #|| $base eq 'N' 
            ) {
            return 0;
        }
    }
    return 1;
}

sub variant_for_base {
    my ($class, $position_info, $base) = @_;
    if($position_info->{variant} eq $base) {
        return $position_info;
    }
    return;
}

sub variants_in_read {
    my ($class, $variant_array, @bases) = @_;
    my @variant_vector = grep { $_ } pairwise { $class->variant_for_base($a, $b) } @$variant_array, @bases;
    return @variant_vector;
}
        
