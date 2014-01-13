#!/usr/bin/perl -w
use strict;

my @result = `qstat | grep y535`;

while(@result){
	$_ = pop @result;
	/(\d+)\.iox/;
	print  "qdel $1\n";
	system "qdel $1";
}
