default_target: compile
	

build:
	mkdir -p build

build/Makefile: build
	cd build && cmake ..

build/src/main: build/Makefile
	cd build && make


# available commands:
compile: build/src/main

run: build/src/main
	./build/src/main

clean:
	cd build && make clean

cleanall:
	rm -rf build/*

help:
	@echo "Hi, this is FPSR main application build system."
	@echo
	@echo "List of available commands:"
	@echo " - compile: builds the appliaction in the <build> directory,"
	@echo "            creates the directory if necessary"
	@echo " - run: runs the application, compiles it if needed"
	@echo " - clean: removes the compilation products"
	@echo " - cleanall: ereases build directory"
