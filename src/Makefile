CC		= g++-mp-4.7
GSL		= /opt/local/lib/libgsl.a /opt/local/lib/libgslcblas.a
CC_FLAGS	= -Wall -O3 -DHAVE_INLINE -static-libgcc -static-libstdc++ -fopenmp
LD_FLAGS	= ${GSL} -fopenmp
filterHDobj	= emission.o common-functions.o jump-diffusion.o minimization.o
cloneobj	= clone.o clone-bulk.o clone-prior.o clone-llh.o clone-predict.o clone-update.o clone-fwd-bwd.o
cloneHDobj	= $(cloneobj) emission.o common-functions.o minimization.o log-space.o cloneHD-functions.o cloneHD-inference.o
.PHONY: clean all
all: filterHD cloneHD
	rm -f ./*.o
filterHD: $(filterHDobj) filterHD.o
	$(CC) $(CC_FLAGS) $^ -o ../build/filterHD $(LD_FLAGS)
cloneHD: $(cloneHDobj) cloneHD.o
	$(CC) $(CC_FLAGS) $^ -o ../build/cloneHD $(LD_FLAGS)	
%.o: %.cpp
	$(CC) $(CC_FLAGS) -c $< -o $@
$(cloneobj): clone.h $(cloneobj:.o=.cpp)
$(filterHDobj) $(cloneHDobj): $(.TARGET:.o=.h)
clean:
	rm -f ./*.o
