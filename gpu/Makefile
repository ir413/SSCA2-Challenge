include Makefile.var

TARGET = SSCA2

OBJS = main.o

all: sprng $(OBJS)
	$(LD) $(LDFLAGS) $(GENCODE_FLAGS) -o $(TARGET) $(OBJS) $(LIB)

main.o: main.cu
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

sprng:
	(cd sprng2.0; $(MAKE); cd ..)

cleansprng:
	(cd sprng2.0; $(MAKE) clean; cd ..)

clean:  cleansprng
	rm -f *.o *~ $(TARGET)
