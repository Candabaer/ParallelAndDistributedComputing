all: main

main.o: main.cpp
	c++ -std=c++11 -pthread -c $<
main: main.o
	c++ -std=c++11 -o main main.o -pthread

