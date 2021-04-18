all: exec

exec: build/exercise1
	./build/exercise1

build/exercise1: exercise1.cpp | build
	g++ -O3 -ffast-math -march=native exercise1.cpp -o build/exercise1

build: 
	mkdir -p build

clean:
	rm -rf build

.PHONY: all exec clean