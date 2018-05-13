CPPC=g++
CFLAGS=-Wall -Werror -g -O0 -fno-inline -std=c++11

%.o: %.cpp
	$(CPPC) -c -o $@ $< $(CFLAGS)

clean:
	rm -rf *.o relaxlag resultValidator resultAggregator

relaxlag: relaxlag.o
	$(CPPC) -o $@ $^ $(CFLAGS)

resultValidator: resultValidator.o
	$(CPPC) -o $@ $^ $(CFLAGS)

resultAggregator: resultAggregator.o
	$(CPPC) -o $@ $^ $(CFLAGS)

