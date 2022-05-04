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
testmap: maptest
	./maptest 100000 100
sorttest: sorttest.cpp samplesortlib.hpp
	g++ $(CPPFLAGS) $(OMPFLAGS)  -DAVX512 -o sorttest sorttest.cpp
sorttestcpp17: sorttestcpp17.cpp samplesortlib.hpp
	g++ $(CPPFLAGS) $(OMPFLAGS)  -DAVX512 -o sorttestcpp17 sorttestcpp17.cpp
sorttest.single : sorttest.cpp samplesortlib.hpp
	g++ $(CPPFLAGS)  -o sorttest.single  sorttest.cpp
maptest: maptest.cpp  samplesortlib.hpp simplemap.hpp
	g++ $(CPPFLAGS) $(OMPFLAGS)  -o maptest maptest.cpp
sorttest128.a64fx : sorttest128.cpp samplesortlib.hpp
	FCC    -g -fopenmp -DSVE -Nclang  -Ofast -mcpu=a64fx+sve -DSORTLIB_MEASURE_TIME  -o sorttest128.a64fx  sorttest128.cpp 
sorttest.a64fx : sorttest.cpp samplesortlib.hpp
	FCC    -g -fopenmp -DSVE -Nclang  -Ofast -mcpu=a64fx+sve -DSORTLIB_MEASURE_TIME -I ../simdsort -o sorttest.a64fx  sorttest.cpp 
sorttest128: sorttest128.cpp samplesortlib.hpp
	g++  -std=c++1z $(CPPFLAGS) -DSORTLIB_MEASURE_TIME -DAVX512   -o sorttest128  sorttest128.cpp 
