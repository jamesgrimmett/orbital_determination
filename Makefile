CC = g++

# compiler flags:
CFLAGS = -g -Wall -std=c++11

all: lls_ex10_1 lls_ex10_2

lls_ex10_1: ex10-1_linear_least_squares.cc
	$(CC) $(CFLAGS) -o lls_1 ex10-1_linear_least_squares.cc

lls_ex10_2: ex10-2_linear_least_squares.cc
	$(CC) $(CFLAGS) -o lls_2 ex10-2_linear_least_squares.cc
