#!/usr/bin/env perl

use Test::Most;
use FindBin;
plan tests => 4;
bail_on_fail if 0;
use Env::Path 'PATH';

if (-d $FindBin::RealBin) {
  $ENV{PATH} .= ":".$FindBin::RealBin;
}
ok(scalar PATH->Whence($_), "$_ in PATH") for qw(blastn bedtools kraken parallel);
