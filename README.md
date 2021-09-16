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

Class T should provide a member function getsortkey()
using which this function sort the array.


## Limitations

Can fail if the size of the array is smaller than the square of the
number of threads.

Assumes that the stack size is large enough to place working arrays
(n*(sizeof(T)+32) bytes). Use ulimit -h (or limit stacksize in Csh),
or modify the code so that it uses new/delete. 

## Performance sample

![Performance on Fugaku. N is the number of elements, Sort time is in second](fugaku.jpg)

Measured on Fugaku using sorttest.cpp. N is the number of
elements. Sort time is in second. The dashed curve with filled
triangles is the time for std::sort. Filled squares, pentagons,
open triangles, squares and pentagons are the results of
samplesort_bodies called with 2, 4,  12,
24, 48 threads.

## Algorithm

Sample sort with sorting applied to index-value pairs (locally
generated). The original array is reordered according to the sorted
index-value pair. Thus this library is optimized for classes with
relatively large sizes (more than 32 bytes).
