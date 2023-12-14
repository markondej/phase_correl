FLAGS = -Wall -O3
LIBS = -lm

all:
	g++ $(FLAGS) -std=c++11 -o phase_correl_cpp phase_correl.cpp $(LIBS)
	gcc $(FLAGS) -o phase_correl_c phase_correl.c $(LIBS)

