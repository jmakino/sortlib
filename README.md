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


