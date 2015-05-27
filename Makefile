# Makefile for simulation

# The following flags work for OSX.
LDFLAGS = -framework GLUT -framework OpenGL

# For linux, you must link to the GLUT and OpenGL libraries
# and use the following flags
#LDFLAGS=-lglut -lGLU -lGL -lGLEW 

# Suppresses warnings about the deprecated GLUT routines.
CFLAGS = -Wno-deprecated -O3 -march=native

all: galaxy_simulator barnes_hut

barnes_hut: barnes_hut.o graphics.o
	gcc -o barnes_hut barnes_hut.o graphics.o $(LDFLAGS)

barnes_hut.o: barnes_hut.c
	gcc -c barnes_hut.c $(CFLAGS)

galaxy_simulator: galaxy_simulator.o graphics.o
	gcc -o galaxy_simulator galaxy_simulator.o graphics.o $(LDFLAGS)

galaxy_simulator.o: galaxy_simulator.c
	gcc -c galaxy_simulator.c $(CFLAGS)
	
graphics.o: graphics.c
	gcc -c graphics.c $(CFLAGS)	
	
clean:
	rm -f ./galaxy_simulator ./barnes_hut *.o
