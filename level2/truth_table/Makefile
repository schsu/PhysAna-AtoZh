SRC := $(wildcard *.cpp)
OBJ := $(SRC:.cpp=.o)
CFLAGS := -O3
LDFLAGS= -L$(DELPHESPATH) -lDelphes
INCLUDES = -I$(DELPHESPATH)
ROOTCFLAGS := $(shell root-config --cflags)
ROOTGLIBS := $(shell root-config --glibs)

%.o:%.cpp
	g++ -c -o $@ $< $(CFLAGS) $(INCLUDES) $(ROOTCFLAGS)

truth_table: $(OBJ)
	g++ $(OBJ) $(CFLAGS) $(LDFLAGS) $(ROOTGLIBS) -o truth_table

clean:
	rm -rf *.o truth_table
