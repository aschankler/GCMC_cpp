ROOT_DIR=$(shell pwd)
ODIR  = $(ROOT_DIR)/obj
SDIR  = $(ROOT_DIR)/src

CXX   = g++
CFLAG = -std=c++11
 
DEPS  = $(shell ls $(SDIR)/*.h)
SRC   = $(shell ls $(SDIR)/*.cpp)
OBJ   = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRC))

gcmc.x : $(OBJ)
	$(CXX) -o $@ $^ $(CFLAG)

$(ODIR)/%.o : $(SDIR)/%.cpp $(DEPS) | $(ODIR)/.
	$(CXX) -c -o $@ $< $(CFLAG)

%/. : 
	mkdir -p $(patsubst %/.,%,$@)
	
.PRECIOUS: %/.
.PHONY: clean clean_test clean_all
clean:
	rm -rf gcmc.x $(ODIR)
clean_test:
	rm -rf *.axsf qe.in qe.out log.dat
clean_all:
	rm -rf gcmc.x $(ODIR) *.axsf qe.in qe.out log.dat
