CC = g++
CXXFLAGS = -std=c++11 -pg

program     = binary_disruption
objects     = binary.o encounters.o MC_velocity.o random_direction.o random_numbers.o vector_maths.o main.o

$(program): $(objects)
	$(CC) $(objects) -o $(program) $(CXXFLAGS)
	
binary.o: binary.cpp constants.h vector_maths.h random_numbers.h
	$(CC) -c binary.cpp $(CXXFLAGS)

encounters.o: encounters.cpp constants.h vector_maths.h random_numbers.h random_direction.h binary.h MC_velocity.h
	$(CC) -c encounters.cpp $(CXXFLAGS)

MC_velocity.o: MC_velocity.cpp random_numbers.h
	$(CC) -c MC_velocity.cpp $(CXXFLAGS)

random_direction.o: random_direction.cpp constants.h vector_maths.h random_numbers.h
	$(CC) -c random_direction.cpp $(CXXFLAGS)

random_numbers.o: random_numbers.cpp
	$(CC) -c random_numbers.cpp $(CXXFLAGS)

vector_maths.o: vector_maths.cpp
	$(CC) -c vector_maths.cpp $(CXXFLAGS)

main.o: main.cpp constants.h random_direction.h random_numbers.h vector_maths.h MC_velocity.h binary.h encounters.h
	$(CC) -c main.cpp $(CXXFLAGS)

clean:
	rm -f  core $(objects) $(program)
