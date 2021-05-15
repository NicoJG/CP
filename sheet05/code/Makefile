all: exercise2 exercise3 

exercise2: build/exercise2
	./build/exercise2
	python exercise2.py

build/exercise2: exercise2.cpp | build
	g++ -O3 -ffast-math -march=native exercise2.cpp -o build/exercise2

exercise3: build/exercise3
	./build/exercise3
	python exercise3.py

build/exercise3: exercise3.cpp | build
	g++  -O3 -ffast-math -march=native exercise3.cpp -o build/exercise3

build: 
	mkdir -p build

clean:
	rm -rf build

.PHONY: all  exercise2 exercise3 clean 