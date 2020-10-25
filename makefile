FLAGS = -Wall -O3
LIBS = -lm

all:
	g++ $(FLAGS) $(LIBS) -o phase_correl phase_correl.cpp

clean:
	rm ./phase_correl
