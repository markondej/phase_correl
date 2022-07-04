FLAGS = -Wall -O3
LIBS = -lm

all:
	g++ $(FLAGS) -o phase_correl_cpp phase_correl.cpp $(LIBS)
	gcc $(FLAGS) -o phase_correl_c phase_correl.c $(LIBS)

clean:
	rm ./phase_correl
