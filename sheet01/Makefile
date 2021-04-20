all: code latex

latex: 
	make -C latex

code:
	make -C code

clean:
	make -C latex clean
	make -C code clean

.PHONY: all latex code clean