CPPC=g++
CFLAGS=-Wall -Werror -g -O0 -fno-inline
DEPS = relaxlag.cpp
OBJ = relaxlag.o

%.o: %.cpp $(DEPS)
	$(CPPC) -c -o $@ $< $(CFLAGS)

relaxlag: $(OBJ)
	$(CPPC) -o $@ $^ $(CFLAGS)

clean:
	rm -rf $(OBJ) relaxlag

