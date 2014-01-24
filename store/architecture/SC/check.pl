while(1){
	system "date";
	system "qstat | grep y535 | head -1";
	system "ls result | tail -1";
	sleep(1);
}


