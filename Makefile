
all: laneemden pressure convergence

laneemden: laneemden.cpp csvwriter
	clang++ laneemden.cpp build/csvwriter.o -O3 -lpthread -std=c++17 -o laneemden

pressure: pressure.cpp csvwriter build/data.o build/polyfit.o
	clang++ pressure.cpp build/csvwriter.o build/data.o build/polyfit.o -O3 -lpthread -std=c++17 -o pressure

convergence: csvwriter
	clang++ convergence.cpp build/csvwriter.o -O3 -lpthread -std=c++17 -o error

csvwriter: csvwriter.cpp csvwriter.h
	clang++ -c csvwriter.cpp -O3 -std=c++17 -o build/csvwriter.o

build/data.o: data.cpp data.h
	clang++ -c data.cpp -O3 -std=c++17 -o build/data.o

build/polyfit.o: polyfit.cpp polyfit.h
	clang++ -c polyfit.cpp -O3 -std=c++17 -o build/polyfit.o

clean:
	rm -f ./laneemden build/* out/* ./error ./pressure

dirs:
	mkdir build out
