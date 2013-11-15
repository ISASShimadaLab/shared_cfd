while(1){
	system "date";
	system "qstat | grep y535";
	system "ls result | tail -1";
	sleep(1);
}


