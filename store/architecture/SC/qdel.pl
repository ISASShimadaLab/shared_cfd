#!/usr/bin/perl -w
use strict;

my @result = `qstat | grep y535`;

while(@result){
	$_ = pop @result;
	/(\d+)\.jxclstm/;
	print  "qdel -k $1\n";
	system "qdel -k $1";
}
