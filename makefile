BINS = diffusion 
LIBRARIES = -L/usr/lib64/openmpi/lib -lm -lmpi -Wl,-rpath,/usr/lib64/openmpi/lib -Wl,--enable-new-dtags
CFLAGS  = -O3 -I. -I/usr/include/openmpi-x86_64 -pthread -fopenmp $(LIBRARIES)

.PHONY : all
all : $(BINS) 

diffusion: diffusion.c
	gcc $(CFLAGS) -o $@ $<

.PHONY : clean
clean :
	rm -rf $(BINS)
