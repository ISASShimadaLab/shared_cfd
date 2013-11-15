while(1){
	system "date";
	my @statement =`qstat | grep y535`;
        if(@statement != 0){
                print @statement[0];
        }else{
                open PP,"| sh mpirun.sh";
                print PP "\n";
                close PP;
        }
	system "ls result | tail -1";
	sleep(1);
}


