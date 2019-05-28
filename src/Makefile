CXX   = g++
CFLAG = -std=c++11
ODIR  = ../obj
 
DEPS  = $(shell ls *.h)
SRC   = $(shell ls *.cpp)
OBJ   = $(patsubst %.cpp,$(ODIR)/%.o,$(SRC))

$(shell mkdir -p $(ODIR))

../gcmc.x : $(OBJ)
	$(CXX) -o $@ $^ $(CFLAG)

$(ODIR)/%.o : %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAG)
	
.PHONY: clean

clean:
	rm -rf ../gcmc.x $(ODIR)
