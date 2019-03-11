
all: laneemden1 pressure convergence

laneemden1: laneemden1.cpp csvwriter
	clang++ laneemden1.cpp build/csvwriter.o -O3 -lpthread -std=c++17 -o laneemden1

laneemden2: laneemden2.cpp csvwriter
	clang++ laneemden2.cpp build/csvwriter.o -O3 -lpthread -std=c++17 -o laneemden2

pressure: pressure.cpp csvwriter build/data.o build/polyfit.o
	clang++ pressure.cpp build/csvwriter.o build/data.o build/polyfit.o -O3 -lpthread -std=c++17 -o pressure

convergence:
	clang++ convergence.cpp -O3 -lpthread -std=c++17 -o error

csvwriter: csvwriter.cpp csvwriter.h
	clang++ -c csvwriter.cpp -O3 -std=c++17 -o build/csvwriter.o

build/data.o: data.cpp data.h
	clang++ -c data.cpp -O3 -std=c++17 -o build/data.o

build/polyfit.o: polyfit.cpp polyfit.h
	clang++ -c polyfit.cpp -O3 -std=c++17 -o build/polyfit.o


clean:
	rm -f ./laneemden1 ./laneemden2 build/* out/* ./error ./pressure

dirs:
	mkdir build out
