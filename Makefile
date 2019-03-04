
all: laneemden1 laneemden2 pressure

laneemden1: laneemden1.cpp csvwriter
	clang++ laneemden1.cpp build/csvwriter.o -O3 -lpthread -std=c++17 -o laneemden1

laneemden2: laneemden2.cpp csvwriter
	clang++ laneemden2.cpp build/csvwriter.o -O3 -lpthread -std=c++17 -o laneemden2

pressure: pressure.cpp csvwriter
	clang++ pressure.cpp build/csvwriter.o -lpthread -std=c++17 -o pressure

csvwriter: csvwriter.cpp csvwriter.h
	clang++ -c csvwriter.cpp -O3 -std=c++17 -o build/csvwriter.o

clean:
	rm -f ./laneemden1 ./laneemden2 build/*

