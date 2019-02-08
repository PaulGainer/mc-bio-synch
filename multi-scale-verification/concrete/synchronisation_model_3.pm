dtmc

// Sensor states
const int protocol_start		= -1;
const int start 			= 0;
const int update			= 1;
const int maxModeValue			= update;
const int T 			        = 10;

const int RP;
const double mu;
const double epsilon ;

const N = 3;

global count:[0..N] ;



init 
	(count = 0) &
	(s1Mode = protocol_start) & 
	(s2Mode = protocol_start) & 
	(s3Mode = protocol_start) 
endinit


formula pert1 = (s1Phase*epsilon*count - floor(s1Phase * epsilon *count)) >= 0.5 ? ceil(s1Phase*epsilon*count) : floor(s1Phase*epsilon*count);
formula pert2 = (s2Phase*epsilon*count - floor(s2Phase * epsilon *count)) >= 0.5 ? ceil(s2Phase*epsilon*count) : floor(s2Phase*epsilon*count);
formula pert3 = (s3Phase*epsilon*count - floor(s3Phase * epsilon *count)) >= 0.5 ? ceil(s3Phase*epsilon*count) : floor(s3Phase*epsilon*count);

module env 
	envMode :[0..maxModeValue] ;
	[] (envMode = start) -> (count' = 0) & (envMode' = update);
	[updateClock] (envMode = update) -> (envMode' = start);
endmodule


module sensor1	
	s1Mode : [-1..maxModeValue] ;
 	s1Phase : [1..T] ;
	[protocol_start] (s1Mode = protocol_start) -> (s1Mode' = start);
// START
	[] 	  (s1Mode = start) 
		& (s1Phase = T) 
		& (envMode = update)
	    	& (count < N)
	-> 
		(1-mu): (s1Mode' = update) & (count' = count+1) 
		+(mu) : (s1Mode' = update) ;
	[]   	  (s1Mode = start) 
	  	& (s1Phase + pert1+1 > T)
		& (envMode = update)
		& (count < N)
		& (s1Phase != T)
		& (s2Phase = T | s3Phase = T)
		& (s2Mode = update | s2Phase <= s1Phase)
		& (s3Mode = update | s3Phase <= s1Phase)
	-> 
		(1-mu): (s1Mode' = update) & (count' = count+1) 
		+(mu) : (s1Mode' = update) ;

	[]   	  (s1Mode = start) 
	  	& (s1Phase + pert1+1  <= T)
		& (envMode = update)
		& (count < N)
		& (s1Phase != T)
		& (s2Phase = T | s3Phase = T)
		& (s2Mode = update | s2Phase <= s1Phase)
		& (s3Mode = update | s3Phase <= s1Phase)
	->
		(s1Mode' = update);
	[]   	  (s1Mode = start) 
		& (s1Phase != T) 
		& (s2Phase != T)
		& (s3Phase != T)
		& (envMode=update)  
	-> 							
		(s1Mode' = update) ;
// UPDATE CLOCK & SYNC
	[sync]	  (s1Mode = update) 
			& (s1Phase = T) 						
	-> 
			  (s1Mode' = start) 
			& (s1Phase' = 1);

	[sync] 	  (s1Mode = update) 
			& (s1Phase < T) 
		    	& (s1Phase <= RP)							
	-> 
		  	  (s1Mode' = start) 
			& (s1Phase' = s1Phase + 1) ;
	[sync]     (s1Mode = update) 
			& (s1Phase < T)  
		       	& (s1Phase > RP)	
			& (s1Phase + pert1 + 1 <= T)	
	-> 
			  (s1Mode' = start) 
			& (s1Phase' = s1Phase + pert1 +1 );
	[sync]	  (s1Mode = update) 
			& (s1Phase < T) 
			& (s1Phase > RP)	
			& (s1Phase + pert1 + 1 > T)	
	-> 
	   		  (s1Mode' = start) 
			& (s1Phase' = 1);
endmodule

module sensor2	
	s2Mode : [-1..maxModeValue] ;
 	s2Phase : [1..T] ;
	[protocol_start] (s2Mode = protocol_start) -> (s2Mode' = start);

	// START
	[]	  (s2Mode = start) 
		& (s2Phase = T) 
		& (envMode = update)
		& (count < N)	
	-> 
		(1-mu): (s2Mode' = update) & (count' = count+1) 
		+(mu) : (s2Mode' = update) ;

	[]   	  (s2Mode = start) 
	  	& (s2Phase + pert2 +1> T)
		& (envMode = update)
		& (count < N)
		& (s2Phase != T)
		& (s1Phase = T | s3Phase = T)
		& (s1Mode = update | s1Phase <= s2Phase)
		& (s3Mode = update | s3Phase <= s2Phase)
	-> 
		(1-mu): (s2Mode' = update) & (count' = count+1) 
		+(mu) : (s2Mode' = update) ;

	[]   	  (s2Mode = start) 
	  	& (s2Phase + pert2 +1 <= T)
		& (envMode = update)
		& (count < N)
		& (s2Phase != T)
		& (s1Phase = T | s3Phase = T)
		& (s1Mode = update | s1Phase <= s2Phase)
		& (s3Mode = update | s3Phase <= s2Phase)
	-> 
		(s2Mode' = update)
	;
	[] 	  (s2Mode = start) 
		& (s1Phase != T) 
		& (s2Phase != T)
		& (s3Phase != T)
		& (envMode=update)  
	-> 
	 	  (s2Mode' = update) ;

	// UPDATE CLOCK & SYNC
	[sync] (s2Mode = update) & (s2Phase = T) 						-> 
							(s2Mode' = start) & (s2Phase' = 1);
	
	[sync] 	  (s2Mode = update) & (s2Phase < T) 
		    	& (s2Phase <= RP)							
		-> 
		  	  (s2Mode' = start) 
			& (s2Phase' = s2Phase + 1) 
			;
	[sync]    (s2Mode = update) & (s2Phase < T)  
		       	& (s2Phase > RP)	
			& (s2Phase + pert2 +1 <= T)
		-> 
			  (s2Mode' = start) 
			& (s2Phase' = s2Phase + pert2 +1)
							;
	[sync]	  (s2Mode = update) 
			& (s2Phase < T) 
			& (s2Phase > RP)	
			& (s2Phase +pert2 + 1 > T)	
		-> 
	   		  (s2Mode' = start) & (s2Phase' = 1) 				
			;

endmodule

module sensor3	
	s3Mode : [-1..maxModeValue] ;
 	s3Phase : [1..T] ;
	[protocol_start] (s3Mode = protocol_start) -> (s3Mode' = start);

	// START
	[] 	  (s3Mode = start) 
		& (s3Phase = T) 
		& (envMode = update)
		& (count < N)	
	-> 
							(1-mu): (s3Mode' = update) & (count' = count+1) 
							+(mu) : (s3Mode' = update) ;

	[]   	  (s3Mode = start) 
	  	& (s3Phase + pert3 +1  > T)
		& (envMode = update)
		& (count < N)
		& (s3Phase != T)
		& (s1Phase = T | s2Phase = T)
		& (s1Mode = update | s1Phase <= s3Phase)
		& (s2Mode = update | s2Phase <= s3Phase)
	-> 
		(1-mu): (s3Mode' = update) & (count' = count+1) 
		+(mu) : (s3Mode' = update) ;

	[]   	  (s3Mode = start) 
	  	& (s3Phase + pert3  +1 <= T)
		& (envMode = update)
		& (count < N)
		& (s3Phase != T)
		& (s1Phase = T | s2Phase = T)
		& (s1Mode = update | s1Phase <= s3Phase)
		& (s2Mode = update | s2Phase <= s3Phase)
	-> 
		  (s3Mode' = update) ;	

	[] (s3Mode = start)  
		& (s1Phase != T) 
		& (s2Phase != T)
		& (s3Phase != T)
		& (envMode=update)  
-> 
							(s3Mode' = update) ;

// UPDATE CLOCK & SYNC
	[sync] (s3Mode = update) & (s3Phase = T) 						-> 
							(s3Mode' = start) & (s3Phase' = 1);
	
	[sync] 	  (s3Mode = update) & (s3Phase < T) 
		    	& (s3Phase <= RP)							
		-> 
		  	  (s3Mode' = start) 
			& (s3Phase' = s3Phase + 1) 
			;
	[sync]    (s3Mode = update) & (s3Phase < T)  
		       	& (s3Phase > RP)	
			& (s3Phase + pert3 + 1 <= T)	
		-> 
			  (s3Mode' = start) 
			& (s3Phase' = s3Phase + pert3 + 1) 
							;
	[sync]	  (s3Mode = update) 
			& (s3Phase < T) 
			& (s3Phase > RP)	
			& (s3Phase +  pert3 + 1 > T)
		-> 
	   		  (s3Mode' = start) & (s3Phase' = 1) 				
			;
endmodule

