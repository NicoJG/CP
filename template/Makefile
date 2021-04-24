all: code latex

code:
	make -C code

latex: 
	make -C latex

clean:
	make -C latex clean
	make -C code clean

.PHONY: all latex code clean