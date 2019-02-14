
all: laneemden1 laneemden2

laneemden1: laneemden1.cpp
	clang++ laneemden1.cpp -std=c++17 -o laneemden1

laneemden2: laneemden2.cpp
	clang++ laneemden2.cpp -std=c++17 -o laneemden2

clean:
	rm -f ./laneemden1 ./laneemden2

