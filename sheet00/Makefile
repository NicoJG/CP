all: exec

exec: build/ex2 build/ex3
	./build/ex2
	python ex2.py
	./build/ex3
	python ex3.py

build/ex2: ex2.cpp | build
	g++ -O3 -ffast-math -march=native ex2.cpp -o build/ex2

build/ex3: ex3.cpp | build
	g++ -O3 -ffast-math -march=native ex3.cpp -o build/ex3

build: 
	mkdir -p build

clean: 
	rm -rf build

.PHONY: all exec clean