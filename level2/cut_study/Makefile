SRC := $(wildcard *.cpp)
OBJ := $(SRC:.cpp=.o)
CFLAGS := -O3
LDFLAGS= -L$(DELPHESPATH) -lDelphes
INCLUDES =  -I$(DELPHESPATH)
ROOTCFLAGS := $(shell root-config --cflags)
ROOTGLIBS := $(shell root-config --glibs)

%.o:%.cpp
	g++ -c -o $@ $< $(CFLAGS) $(INCLUDES) $(ROOTCFLAGS)

cut_study: $(OBJ)
	g++ $(OBJ) $(CFLAGS) $(LDFLAGS) $(ROOTGLIBS) -o cut_study

clean:
	rm -rf *.o cut_study
