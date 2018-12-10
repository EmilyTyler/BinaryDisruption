LIBS= -I/usr/local/include/ -L/usr/lib64 -lgsl -lgslcblas -lm -std=c++0x -O3
CCOMP = g++
LINKER = -o 
CFLAGS = -c

program     = binary_disruption
objects     = binary.o encounters.o MC_velocity.o random_direction.o random_numbers.o vector_maths.o main.o

$(program): $(objects)
	g++ $(objects) -o $(program) $(LIBS)
	
binary.o: binary.cpp constants.h vector_maths.h random_numbers.h
        g++ -c binary.cpp $(LIBS)

encounters.o: encounters.cpp constants.h vector_maths.h random_numbers.h random_direction.h binary.h MC_velocity.h
        g++ -c encounters.cpp $(LIBS)

MC_velocity.o: MC_velocity.cpp random_numbers.h
        g++ -c MC_velocity.cpp $(LIBS)

random_direction.o: random_direction.cpp constants.h vector_maths.h random_numbers.h
        g++ -c random_direction.cpp $(LIBS)

random_numbers.o: random_numbers.cpp
        g++ -c random_numbers.cpp $(LIBS)

vector_maths.o: vector_maths.cpp
        g++ -c vector_maths.cpp $(LIBS)
	
main.o: main.cpp random_direction.h
	g++ -c main.cpp $(LIBS)


clean:
	rm -f  core $(objects) $(program)

#the make file is of the form
#   target: source
#   command
#
#this means that target depends on source, and if source has changed
#then command is implemented, creating a new target. 
#command must be preceeded by a tab.

