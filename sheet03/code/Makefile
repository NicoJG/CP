all: exercise1 exercise2

exercise1: build/exercise1
	./build/exercise1
	python exercise1.py

build/exercise1: exercise1.cpp | build
	g++ -O3 -ffast-math -march=native exercise1.cpp -o build/exercise1

exercise2: build/exercise2
	./build/exercise2
	python exercise2.py

build/exercise2: exercise2.cpp | build
	g++ exercise2.cpp -o build/exercise2

build: 
	mkdir -p build

clean:
	rm -rf build

.PHONY: all exercise1 exercise2 clean