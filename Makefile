% : %.cpp
	g++ $(CPPFLAGS) -o $@ $< 
#CPPFLAGS =  -g   -fsanitize=address # -fsanitize=leak
CPPFLAGS += -fopenmp  -O3 -march=native
#CPPFLAGS += -ffast-math    -ftree-vectorize  -fopt-info-vec-optimized=vector.txt
testsort: sorttest
	echo Test ssamplesort for n=2 to 8000. This will take a while...
	./sorttest 2 8000
	./sorttest 8001
sorttest: sorttest.cpp samplesortlib.hpp
