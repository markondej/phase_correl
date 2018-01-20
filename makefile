FLAGS = -Wall -O3
LIBS = -lm

all:
	gcc $(FLAGS) $(LIBS) -o phase_correl main.c

clean:
	rm ./phase_correl
