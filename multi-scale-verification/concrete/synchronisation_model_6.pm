dtmc

// Sensor states
const int protocol_start		= -1;
const int start 			= 0;
const int update			= 1;
const int maxModeValue			= update;

const int T;
const int RP;
const double mu;
const double epsilon ;

const N = 6;

global count:[0..N] ;



init 
	(count = 0) &
	(s1Mode = protocol_start) & 
	(s2Mode = protocol_start) & 
	(s3Mode = protocol_start) &
	(s4Mode = protocol_start) &
	(s5Mode = protocol_start) &
	(s6Mode = protocol_start)

endinit


formula pert1 = (s1Phase*epsilon*count - floor(s1Phase * epsilon *count)) >= 0.5 ? ceil(s1Phase*epsilon*count) : floor(s1Phase*epsilon*count);
formula pert2 = (s2Phase*epsilon*count - floor(s2Phase * epsilon *count)) >= 0.5 ? ceil(s2Phase*epsilon*count) : floor(s2Phase*epsilon*count);
formula pert3 = (s3Phase*epsilon*count - floor(s3Phase * epsilon *count)) >= 0.5 ? ceil(s3Phase*epsilon*count) : floor(s3Phase*epsilon*count);
formula pert4 = (s4Phase*epsilon*count - floor(s4Phase * epsilon *count)) >= 0.5 ? ceil(s4Phase*epsilon*count) : floor(s4Phase*epsilon*count);
formula pert5 = (s5Phase*epsilon*count - floor(s5Phase * epsilon *count)) >= 0.5 ? ceil(s5Phase*epsilon*count) : floor(s5Phase*epsilon*count);
formula pert6 = (s6Phase*epsilon*count - floor(s6Phase * epsilon *count)) >= 0.5 ? ceil(s6Phase*epsilon*count) : floor(s6Phase*epsilon*count);

module env
	envMode :[0..maxModeValue] ;
	[] (envMode = start) -> (count' = 0) & (envMode' = update);
	[sync] (envMode = update) -> (envMode' = start);
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
		& (s2Phase = T | s3Phase = T | s4Phase = T | s5Phase = T | s6Phase = T)
		& (s2Mode = update | s2Phase <= s1Phase)
		& (s3Mode = update | s3Phase <= s1Phase)
		& (s4Mode = update | s4Phase <= s1Phase)
		& (s5Mode = update | s5Phase <= s1Phase)
		& (s6Mode = update | s6Phase <= s1Phase)
	-> 
		(1-mu): (s1Mode' = update) & (count' = count+1) 
		+(mu) : (s1Mode' = update) ;

	[]   	  (s1Mode = start) 
	  	& (s1Phase + pert1+1  <= T)
		& (envMode = update)
		& (count < N)
		& (s1Phase != T)
		& (s2Phase = T | s3Phase = T | s4Phase = T | s5Phase = T | s6Phase = T)
		& (s2Mode = update | s2Phase <= s1Phase)
		& (s3Mode = update | s3Phase <= s1Phase)
		& (s4Mode = update | s4Phase <= s1Phase)
		& (s5Mode = update | s5Phase <= s1Phase)
		& (s6Mode = update | s6Phase <= s1Phase)
	->
		(s1Mode' = update);
	[]   	  (s1Mode = start) 
		& (s1Phase != T) 
		& (s2Phase != T)
		& (s3Phase != T)
		& (s4Phase != T)
		& (s5Phase != T)
		& (s6Phase != T)
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



module sensor2	= sensor1 [s1Mode = s2Mode, s2Mode = s1Mode, s1Phase = s2Phase, s2Phase = s1Phase] endmodule
module sensor3	= sensor1 [s1Mode = s3Mode, s3Mode = s1Mode, s1Phase = s3Phase, s3Phase = s1Phase] endmodule
module sensor4	= sensor1 [s1Mode = s4Mode, s4Mode = s1Mode, s1Phase = s4Phase, s4Phase = s1Phase] endmodule
module sensor5	= sensor1 [s1Mode = s5Mode, s5Mode = s1Mode, s1Phase = s5Phase, s5Phase = s1Phase] endmodule
module sensor6	= sensor1 [s1Mode = s6Mode, s6Mode = s1Mode, s1Phase = s6Phase, s6Phase = s1Phase] endmodule
