SRC := $(wildcard *.cpp)
OBJ := $(SRC:.cpp=.o)
CFLAGS := -O3
LDFLAGS= -L$(DELPHESPATH) -lDelphes
INCLUDES = -I$(ROOTSYS)/include -I$(DELPHESPATH)

%.o:%.cpp
	g++ -c -o $@ $< $(CFLAGS) $(INCLUDES)

truthstudy: $(OBJ)
	g++ $(OBJ) $(CFLAGS) $(LDFLAGS) -o truthstudy

clean:
	rm -rf *o
