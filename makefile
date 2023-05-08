CC = g++
CXXFLAGS = -std=c++11

EXECUTABLE     = bin/binary_disruption
SOURCES        = $(wildcard src/*.cpp)
OBJECTS        = $(patsubst %.cpp, %.o, $(SOURCES))

all:	build $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXECUTABLE) $(CXXFLAGS)
    
$(OBJECTS): src/%.o : src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

build:
	@mkdir -p bin

clean:
	rm -rf $(EXECUTABLE) $(OBJECTS) bin
