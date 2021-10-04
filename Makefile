% : %.cpp
	g++ $(CPPFLAGS) -o $@ $< 
#CPPFLAGS =  -g   -fsanitize=address # -fsanitize=leak
CPPFLAGS +=   -O3 -march=native
OMPFLAGS += -fopenmp 
#CPPFLAGS += -ffast-math    -ftree-vectorize  -fopt-info-vec-optimized=vector.txt
testsort: sorttest
	echo Test ssamplesort for n=2 to 8000. This will take a while...
	./sorttest 2 8000
	./sorttest 8001
sorttest: sorttest.cpp samplesortlib.hpp
	g++ $(CPPFLAGS) $(OMPFLAGS) -o sorttest sorttest.cpp
sorttest.single : sorttest.cpp samplesortlib.hpp
	g++ $(CPPFLAGS)  -o sorttest.single  sorttest.cpp
