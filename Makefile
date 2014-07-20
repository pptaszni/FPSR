__Main__: build/Makefile
	cd build && make

build:
	mkdir -p build

build/Makefile: build
	cd build && cmake ..
