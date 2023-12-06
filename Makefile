sources = \
	src/main.cpp \
	src/gl_core_3_3.c
libs = \
	-lGL \
	-lGLU \
	-lglut
inc = \
	-Iinclude
outname = base_freeglut

all:
	g++ -std=c++17 $(sources) $(libs) $(inc) -o $(outname)
clean:
	rm $(outname)
