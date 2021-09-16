# sortlib

An experimental code for test/improve the performance of parallel
sort

you can try this with

```
  make -f Makefile.sorttest sorttest
  env OMP_NUM_THREADS=4 sorttest  120
```
The last line before "success!" message would look like
```
Time to sort = 1.3968e-05 0.000201143 5.41974e-05
```
First number indicates the time for std:sort (single thread),
and second and third both show the time for parallel sample sort,
for first and second calls. Since the overhead of multithread
execution is the largest for the first call, I show the times for
first and second calls.

## Usage as header-only library

```
   #include "PATH/OF/THIS/FILE/samplesortlib.hpp"
```


Provide:
```
   template<class T>
   void samplesort_bodies(T * b,
		       int n)
```
   in namespace SampleSortLib

b is the pointer to array (could be vector?) of class T
n is the number of elements

Class T should provide a menber function getsortkey()
using which this function sort the array.


## Limitations

Can fail if the size of the array is smaller than the squre of the
number of threads.
