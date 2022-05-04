# sortlib

An experimental code for test/improve the performance of parallel
sort and parallel map

you can try this with

```
  make sorttest
  env OMP_NUM_THREADS=4 sorttest  1500
```
The output would look like
```
Time to sort = 0.00029669 0.000591706 0.00015421
Parallel sort success for n=1500 nthreads=4!
Parallel sort success!
```
First number indicates the time for std:sort (single thread),
and second and third both show the times for parallel sample sort,
for the first and second calls. Since the overhead of multithread
execution is the largest for the first call, I show the times for
first and second calls.

This result is ontained  on Intel® Core™ i7-1065G7 with gcc version
7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04).  



## Usage as header-only library

```
   #include "simdsortlib.hpp"
   #include "samplesortlib.hpp"
```

with appropriate -I option

Provide:
```
   template<class T>
   void samplesort(T * b, int n)
```
and
```
   template<class T, class GetKey>
   void samplesort(T * b, int n, GetKey getkey)
```
   in namespace SampleSortLib

b is the pointer to array (could be vector?) of class T and
n is the number of elements

To use the first form, class T should provide a member function getsortkey()
using which this function sort the array. In the second form the
function to get the sort key is given as the third argument. For
example, if we want to sort n elements of class B in stored in array a
using the member variable x as key, we write
```
    SampleSortLib::samplesort(a, n,  [](B & b) ->auto{return b.x;} );
```

The old name samplesort_bodies is kept for compatibility.

Samplesort uses std::sort as the core sorter for single-core, in-cache sort. It
is also possible to use SIMDsort  AVX512 or SVE
(https://github.com/jmakino/simdsort), if the GetKey type is
effectively an 128-bit unsigned integer type and support two member
functions:
```
uint64_t get_hi_key();
uint64_t get_lo_key();
```
(see example in sorttest128.cpp). SIMD version requires either
-DAVX512 or -DSVE specified. Speedup might be marginal with AVX512 but
reasonable with SVE (not large, though, around 20-30%)


## Limitations


Assumes that the stack size is large enough to place working arrays
(n*(sizeof(T)+32) bytes). Use ulimit -s (or limit stacksize in Csh),
or modify the code so that it uses new/delete or other memory
management schemes.

Switches from single-thread std::sort to parallel sort at
size=1000. This might not be the optimum position to switch depending
on architecture.

## Performance sample

![Performance on Fugaku. N is the number of elements, Sort time is in second](fugaku.jpg)

Measured on Fugaku using sorttest.cpp. N is the number of
elements. Sort time is in seconds. The dashed curve with filled
triangles is the time for std::sort. Filled squares, pentagons,
open triangles, squares and pentagons are the results of
samplesort_bodies called with 2, 4,  12,
24, 48 threads.

### Comparison with GCC parallel STL

samplesortcpp17.cpp compares the time for std::sort,
SampleSortLib::samplesort and __gnu_parallel::sort. Example
performance(measured with env OMP_NUM_THREADS=16 sorttestcpp17 100000
100004 1 and other runs with different N)

CPU: Intel(R) Xeon(R) Gold 6140 CPU @ 2.30GHz
OS: Cent 7
Compiler: gcc 11.2.1

100K:0.0117821 0.00246628 0.00178517
1M: 0.224061 0.0275078 0.0223829
10M: 2.89478 0.287627 0.487174




## Algorithm

Sample sort with sorting applied to index-value pairs (locally
generated). The original array is reordered according to the sorted
index-value pair. Thus this library is optimized for classes with
relatively large sizes (more than 32 bytes).

# Map replacement

Experimental library to replace std::map with fast paralell algorithm.
This library is not intended as general-purpose replacement of map,
but desgined for a rather specific use case, where we have an array
of values and want to obtain the indices within array from values.


## Usage as header-only library

```
   #include "PATH/OF/THIS/FILE/samplesortlib.hpp"
   #include "PATH/OF/THIS/FILE/simplemap.hpp"
```

## Limitations:

* The size of array n should be specified before setting any values
* One can use values as indices only after all values corresponding to
  all indices in the specified range are set
* For parallel index search, one should call makemap member function before
  use.
* One makemap is called, or index serch is done, one cannot change the
  content of map. One can clear the map and set new values.
* One can use [] operator only for assignment and not for
  reference. For reference, use at() instead.
  
## Supported functions and API

```
namespace SimpleMapLib{
    template <typename KeyT, typename ValT>
    class Map{
        //std::map like APIs
	Map();
	void clear();
	KeyValuePair* find(KeyT k);
	ValT at(KeyT k);
	KeyAddressPair  operator [](KeyT  k);
	KeyValuePair* end();
        //additional APIs
	Map(ValT s);
        void resize(ValT s);
	void makemap();
    };
}
```
One can create a map by, for example,
```
    SimpleMapLib::Map<int64_t, int> m;
```
In this case, one first set the number of element by
```
    m.resize(n);
```
after this, one can set values by
```
m[key]=index
```
but ALL indices of the range [0,n) must be set before the map m is
used to get index from key, and after all indices are set, makemap
function should be called. The setting of values can be done in
parallel, using for example OpenMP parallel for.

To get index from key. one can use at(), but not []. One can also use
find() and first/second members.

The type of return values of functions find() and end() is a pointer
to SimpleMapLib::KeyAddressPair, a class, which can be used in place
of iterator.


## Changelog

### May 4  2022

* Comparison with __gnu_parallel::sort added
 (sorttestcpp17.cpp). Second number gives the sample sort performance
 and the third gun parallel sort. Our sort seems to be faster for
 large N (100K or larger)


    

### Apr 30 2022

* Fast (on some machines) multipartion algorithm used for initial
  division
* Sample sort is also used for single-thread sort. 

    


