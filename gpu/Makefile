include Makefile.var

TARGET = SSCA2
TEST_TARGET = tests

OBJS = main.o \
       BetweennessCentrality.o \
       Configuration.o \
       Graph.o \
       PrefixSums.o \
       ScalableDataGeneration.o \
       Timer.o \
       Validation.o

TEST_OBJS = PrefixSums.o \
	    PrefixSumsTest.o

HEADERS = BetweennessCentrality.h \
          Configuration.h \
	  Graph.h \
	  PrefixSums.h \
	  ScalableDataGeneration.h \
	  Timer.h \
	  Validation.h

TEST_HEADERS = PrefixSums.h

all: Makefile sprng $(OBJS) $(HEADERS)
	$(LD) $(LDFLAGS) $(GENCODE_FLAGS) -o $(TARGET) $(OBJS) $(LIB)

test: $(TEST_OBJS) $(TEST_HEADERS) 
	$(LD) $(LDFLAGS) $(GENCODE_FLAGS) -o $(TEST_TARGET) $(TEST_OBJS) 

main.o: main.cu
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

BetweennessCentrality.o: BetweennessCentrality.cu BetweennessCentrality.h Graph.h
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

Configuration.o: Configuration.cu Configuration.h
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

Graph.o: Graph.cu Graph.h ScalableDataGeneration.h PrefixSums.h
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

PrefixSums.o: PrefixSums.cu PrefixSums.h
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

PrefixSumsTest.o: PrefixSumsTest.cu PrefixSums.h
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

ScalableDataGeneration.o: ScalableDataGeneration.cu ScalableDataGeneration.h Configuration.h
	$(NVCC) $(INC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

Timer.o: Timer.cu Timer.h
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

Validation.o: Validation.cu Validation.h Configuration.h
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

sprng:
	(cd sprng2.0; $(MAKE); cd ..)

cleansprng:
	(cd sprng2.0; $(MAKE) clean; cd ..)

clean:  cleansprng
	rm -f *.o *~ $(TARGET) $(TEST_TARGET)
