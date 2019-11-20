FLAGS = -Wall -O3
LIBS = -lm

all:
	gcc $(FLAGS) $(LIBS) -o phase_correl phase_correl.c

clean:
	rm ./phase_correl
