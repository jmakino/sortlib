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
    auto n= atoi(argv[1]);
    fprintf(stderr, "n=%d\n", n);
    auto bodies = new Body[n];
    auto b2 = new Body[n];
    auto b3 = new Body[n];
    auto b4 = new Body[n];
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
    auto t2a=SampleSortLib::GetWtime();
    SampleSortLib::sort_bodies(b4, n,
			       [](Body & l, Body & r )
			       ->bool{return l.getsortkey()
				       < r.getsortkey();} );
    auto t2=SampleSortLib::GetWtime();
    SampleSortLib::samplesort_bodies(b3, n);
    //				     [](Body & l)
    //				     ->auto{return l.pos[0];} );
    auto t3=SampleSortLib::GetWtime();
    printf("Time to sort = %g %g %g\n", t1-t0, t2a-t1, t3-t2);
    bool ok=true;
    for (auto i=0; i<n; i++){
	if (bodies[i].id != b2[i].id ||
	    bodies[i].id != b4[i].id   ){
	    ok=false;
	    printf("ERROR at:%d  %ld %ldg\n",i,  bodies[i].id, b3[i].id);
	}
    }
    if (ok){
	printf("Parallel sort success!\n");
    }else{
	printf("Parallel sort failed!\n");
    }
	
    exit(0);
    
}

