CPPC=g++
CFLAGS=-Wall -Werror -O3 -fno-inline -std=c++11

%.o: %.cpp
	$(CPPC) -c -o $@ $< $(CFLAGS)

clean:
	rm -rf *.o relaxlag resultValidator resultAggregator outputAggregator

relaxlag: relaxlag.o
	$(CPPC) -o $@ $^ $(CFLAGS)

resultValidator: resultValidator.o
	$(CPPC) -o $@ $^ $(CFLAGS)

resultAggregator: resultAggregator.o
	$(CPPC) -o $@ $^ $(CFLAGS)

outputAggregator: outputAggregator.o
	$(CPPC) -o $@ $^ $(CFLAGS)

