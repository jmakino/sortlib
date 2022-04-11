#include<cmath>
#include<vector>
#include<functional>
#include<algorithm>
#include<exception>
#include<stdexcept>
#include<cassert>
#include<typeinfo>
#include<cstdio>
#include<cstring>
#include<map>
#include<random>
#include<omp.h>
#include<stdlib.h>

#include "simdsort.hpp"

#include "samplesortlib.hpp"


class Body{
public:
    int64_t id;
    double pos[3];
    double vel[3];
    double acc[3];
    double accold[3];
    double pos2[3];
    double vel2[3];
    double phi;
    inline double getsortkey() {return pos[0];}
    
    
};

int body_compare(const void* a, const void*b)
{
    auto sub = ((Body*)a)->getsortkey()-((Body*)b)->getsortkey();
    if (sub < 0.0) {
	return -1;
    }else if (sub > 0.0){
	return 1;
    }else{
	return 0;
    }
}

template<class T>
void sort_bodies_using_qsort(T * b,
		 int n)
{
    qsort(b, n, sizeof(T), body_compare);
}


int main(int argc, char** argv)
{
    std::mt19937 mt;
    //    using namespace SampleSortLib;
    // 一様実数分布
    // [-1.0, 1.0)の値の範囲で、等確率に実数を生成する
    std::uniform_real_distribution<> dist1(-1.0, 1.0);
    auto nstart= atoi(argv[1]);
    auto nend =nstart+1;
    if (argc > 2) nend= atoi(argv[2]);
    auto showtime = false;
    if (argc > 3) {
	showtime = (atoi(argv[3]) == 1);
    }else if (argc <= 2){
	showtime = true;
    }
		    
		    
    
    auto bodies = new Body[nend];
    auto b2 = new Body[nend];
    auto b3 = new Body[nend];
    auto b4 = new Body[nend];
    for (auto n=nstart; n<nend; n++){
	//	printf("n=%d nt=%d\n", n, omp_get_max_threads());
	for (auto i=0; i<n; i++){
	    bodies[i].id=i;
	    for (auto k=0;k<3; k++) bodies[i].pos[k]= dist1(mt);
	}
	// for (auto i=0; i<10; i++){
	// 	printf("pos[%d][0]=%.3g\n", i, bodies[i].pos[0]);
	// }
	for (auto i=0; i<n; i++){
	    b2[i]=bodies[i];
	    b3[i]=bodies[i];
	    b4[i]=bodies[i];
	}
	auto t0=SampleSortLib::GetWtime();
	SampleSortLib::sort_bodies(bodies, n);
	auto t1=SampleSortLib::GetWtime();
	double dummy;
	SampleSortLib::samplesort_bodies(b2, n,
					 [](Body & l)
					 ->auto{return l.pos[0];});
	auto t2=SampleSortLib::GetWtime();
	SampleSortLib::samplesort_bodies(b3, n);
	//				     [](Body & l)
	//				     ->auto{return l.pos[0];} );
	auto t3=SampleSortLib::GetWtime();
	if (showtime) printf("Time to sort = %g %g %g\n", t1-t0, t2-t1, t3-t2);
	bool ok=true;
	for (auto i=0; i<n; i++){
	    if (bodies[i].id != b2[i].id ||
		bodies[i].id != b3[i].id   ){
		ok=false;
		printf("ERROR at:%d  %ld %ldg\n",i,  bodies[i].id, b3[i].id);
		exit(-1);
	    }
	}
	if (ok){
	    if (showtime){
		printf("Parallel sort success for n=%d nthreads=%d!\n", n,
		       SAMPLESORTLIB_OMP_GET_MAX_THREADS);
	    }
	}else{
	    printf("Parallel sort failed!\n");
	}
    }
    printf("Parallel sort success!\n");    
    exit(0);
    
}

