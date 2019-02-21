
all: laneemden1 laneemden2

laneemden1: laneemden1.cpp csvwriter
	clang++ laneemden1.cpp build/csvwriter.o -lpthread -std=c++17 -o laneemden1

laneemden2: laneemden2.cpp csvwriter
	clang++ laneemden2.cpp build/csvwriter.o -lpthread -std=c++17 -o laneemden2

csvwriter: csvwriter.cpp csvwriter.h
	clang++ -c csvwriter.cpp -std=c++17 -o build/csvwriter.o

clean:
	rm -f ./laneemden1 ./laneemden2 build/*

