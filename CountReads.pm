#!/usr/bin/env perl
package CountReads;

__PACKAGE__->run() unless caller();

use strict;
use warnings;
use List::MoreUtils qw(pairwise);
use IO::File;

sub run() {
    my $class = shift;
    my $ifh = shift;
    unless($ifh) {
        print STDERR $ARGV[0], "\n";
        if ($ARGV[0]) {
            $ifh = IO::File->new($ARGV[0], "r")
                or die "Unable to open $ARGV[0]\n";
        }
        else {
            $ifh = IO::File->new_from_fd(fileno(STDIN), "r")
                or die "Unable to open filehandle to STDIN\n";
        }
    }
    my $ofh = shift;
    unless($ofh) {
        $ofh = IO::File->new_from_fd(fileno(STDOUT), "w")
            or die "Unable to open filehandle to STDOUT\n";
    }
    my $variants;
    my %variant_counts;
    my %file_handles;
    my $discarded = 0;
    my $spanning_counts = 0;
    my @position_counts;

    while(<$ifh>) {
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
                $class->print_readname(\%file_handles, $name, $class->correct_readname($read_name));
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
    print $ofh "Variant(s)\tVariant Supporting Reads\tTotal Spanning Reads\tVAF\n";
    for my $variant (sort keys %variant_counts) {
        printf $ofh "%s\t%d\t%d\t%f%%\n", $variant, $variant_counts{$variant}, $spanning_counts, $variant_counts{$variant} / $spanning_counts * 100;
    }
}

sub correct_readname {
    my ($class, $name) = @_;
    $name =~ s/\/ccs$//g;
    return $name;
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

sub create_readname_filename {
    my ($class, $variant_name) = @_;
    my $filename = $variant_name;
    $filename =~ s/,/_/g;
    $filename .= "_readnames.txt";
    return $filename;
}

sub print_readname {
    my ($class, $filehandles_hashref, $variant_name, $readname) = @_;
    my $fh;
    unless(exists $filehandles_hashref->{$variant_name}) {
        my $filename = $class->create_readname_filename($variant_name);
        $fh = IO::File->new($filename,"w") 
            or die "Unable to open $filename to output readnames supporting $variant_name\n";
        $filehandles_hashref->{$variant_name} = $fh;
    }
    else {
        $fh = $filehandles_hashref->{$variant_name};
    }
    print $fh "$readname\n";
}
