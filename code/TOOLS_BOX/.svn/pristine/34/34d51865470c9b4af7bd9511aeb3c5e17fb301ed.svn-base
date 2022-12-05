CC = g++
MCC = g++

LIBS = -lm

all: revolve_binomial revolve_multistage revolve_online no_checkpointing

revolve.o : revolve.cpp 
## $(CC) -c -O3 revolve.cpp
	$(CC) -c -O3 revolve.cpp
use_revolve : revolve.o use_revolve.o 
	$(CC) -o use_revolve revolve.o use_revolve.o $(LIBS) 
use_revolve.o : use_revolve.cpp 

	$(CC) -c -O3 use_revolve.cpp

revolve_binomial : revolve.o revolve_binomial.o 
	$(CC) -o revolve_binomial revolve.o revolve_binomial.o $(LIBS) 
revolve_binomial.o : revolve_binomial.cpp 

	$(CC) -c -O3 revolve_binomial.cpp
	
revolve_multistage : revolve.o revolve_multistage.o 
	$(CC) -o revolve_multistage revolve.o revolve_multistage.o $(LIBS) 
revolve_multistage.o : revolve_multistage.cpp 

	$(CC) -c -O3 revolve_multistage.cpp
	
revolve_online : revolve.o revolve_online.o 
	$(CC) -o revolve_online revolve.o revolve_online.o $(LIBS) 
revolve_online.o : revolve_online.cpp 

	$(CC) -c -O3 revolve_online.cpp
	
no_checkpointing : no_checkpointing.o 
	$(CC) -o no_checkpointing no_checkpointing.o $(LIBS) 
no_checkpointing.o : no_checkpointing.cpp 

	$(CC) -c -O3 no_checkpointing.cpp
## $(CC) -c -O3 use_revolve.cpp
	
clean:
	rm *.o; rm revolve_binomial; rm revolve_multistage;  rm revolve_online;  rm no_checkpointing; 
	rm ./check_rom/*.txt;