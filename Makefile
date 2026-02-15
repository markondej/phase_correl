CROSS_COMPILE =
CC = $(CROSS_COMPILE)gcc
CFLAGS = -Wall -O3
CXX = $(CROSS_COMPILE)g++
CXXFLAGS = $(CFLAGS) -std=c++11
LIBS = -lm

.PHONY: all

all: phase_correl_cpp phase_correl_c

phase_correl_cpp: phase_correl.cpp
	$(CXX) $(CXXFLAGS) -o phase_correl_cpp phase_correl.cpp $(LIBS)

phase_correl_c:
	$(CC) $(CFLAGS) -o phase_correl_c phase_correl.c $(LIBS)
