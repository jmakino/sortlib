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

inline double GetWtime()
{
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC,&ts )){
	printf("GetWtime Failed\n");
    }
    return ((double) ts.tv_sec)+ ts.tv_nsec*1e-9;
}

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

class KeyValuePair{
public:
    int64_t key;
    double value;
    inline double getsortkey() {return value;}    
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
void sort_bodies(T * b,
		 int n)
{
    
    std::sort(b, b+n, 
	      [](T & l, T & r )
	      ->bool{return l.getsortkey() < r.getsortkey();} );
}
template<class T>
void sort_bodies_using_qsort(T * b,
		 int n)
{
    qsort(b, n, sizeof(T), body_compare);
}


void showdt(const char *s,
	      bool initialize=false)
{
#if 1
    static double t0;
    
    if (!initialize){
	printf("%s dt=%g\n",s, GetWtime()-t0);
    }
    t0=GetWtime();
#endif    
}	       


const int randomize_offset = 48271;    
template<class T>
void samplesort_bodies(T * b,
		       int n)
{
    showdt("", true);
    auto nt = omp_get_max_threads();
    //    auto nt = 4;
    auto nwork0= (n+nt-1)/nt;
    auto nsample= nwork0/nt/10;
    if (nsample < 2) nsample = 2;
    auto nsampletotal = nsample * nt;
    double samplearray[nsampletotal];
    int mystart[nt];
    int myend[nt];
    KeyValuePair key[nt][nwork0];
    //    printf("nt=%d  nwork=%d nsample=%d \n", nt, nwork0, nsample);
    int isrcstart[nt][nt+1];
    int ndest[nt][nt];
    int destsize[nt];
    int deststart[nt];
    auto maxdestsize=0;
#pragma omp parallel
    {
	auto it=omp_get_thread_num();
	if (it==0) showdt("Enter parallel region");
	auto mystartdest = it*nsample;
	auto mystartlocal = it*nwork0;
	auto myendlocal =mystartlocal+nwork0;
	if (myendlocal > n) myendlocal=n;
	mystart[it]=mystartlocal;
	myend[it]=myendlocal;
	auto myrange = myendlocal - mystartlocal;
	auto myoffset = randomize_offset ;
	if (myoffset % myrange == 0 ||myrange % myoffset == 0){
	    myoffset ++;
	}
	for(auto ilocal=0; ilocal<nsample; ilocal++){
	    auto myi = ((ilocal + nt)*myoffset)%myrange;
	    //	    printf("it= %d ilocal=%d myi=%d dest=%d source[%d]=%g\n",
	    //		   it, ilocal, myi, mystartdest+ilocal,
	    //		   mystartlocal+myi,b[mystartlocal+myi].getsortkey());
	    
	    samplearray[mystartdest+ilocal]= b[mystartlocal+myi].getsortkey();
	}
	if (it==0) showdt("Sample made");
	//	printf("sort size=%d %d\n", myrange,it);
	for(auto i=0;i<myrange; i++){
	    key[it][i].key =i+mystartlocal;
	    key[it][i].value =b[i+mystartlocal].getsortkey();
	}
	//	sort_bodies(b+mystartlocal, myrange);
	if (it==0) showdt("Key made");
	sort_bodies(key[it], myrange);
	if (it==0) showdt("Initial sort end");
	//	printf("sort end %d\n", it);
	if (it==0){
	    showdt("Intial sampling and internal sort");
	    //    printf("sample sort size=%d\n", nsampletotal);
	    std::sort(samplearray, samplearray+nsampletotal, 
		      [](double  l, double  r )   ->bool{return l < r;} );
	    // for(auto i=0;i<nsampletotal; i+=nsample){
	    // 	printf("sample[%d]=%g\n", i, samplearray[i]);
	    // }
	    showdt(" sample sort");
	}
#pragma omp barrier	
#if 0	
    }
#pragma omp parallel for schedule(static)
    for(auto it=0; it<nt; it++){
	auto mystartdest = it*nsample;
	auto mystartlocal = it*nwork0;
	auto myendlocal =mystartlocal+nwork0;
	if (myendlocal > n) myendlocal=n;
	auto myrange = myendlocal - mystartlocal;
#endif	
	//	auto bstart =  b+mystartlocal;
	auto bstart =  key[it];
	if (myendlocal > n) myendlocal=n;
	isrcstart[it][0]=0;
	isrcstart[it][nt]=myrange;
	// for(auto i=0; i<myrange; i++){
	//     printf(" %g", bstart[i].getsortkey());
	// }
	// printf("\n");
	for (auto itdest=1; itdest<nt; itdest++){
		
	    // find the next start location using bisection
	    double nextval = samplearray[itdest*nsample];
	    int start = isrcstart[it][itdest-1];
	    //	    if (bstart[start].getsortkey() < nextval) start ++;
	    int end = myrange-1;
	    // printf("%d %d %d %d val=%g start,end = %g %g\n",
	    // 	   it, itdest, start, end,
	    // 	   nextval,
	    // 	   bstart[start].getsortkey(),
	    // 	   bstart[end].getsortkey()  );
	    if (isrcstart[it][itdest-1]==myrange){
		//		printf("end already reached\n");
		isrcstart[it][itdest]=myrange;
		
	    }else if (bstart[start].getsortkey() >= nextval){
		//		printf("no element\n");
		isrcstart[it][itdest]=isrcstart[it][itdest-1];
	    }else if (bstart[end].getsortkey() <nextval){
		//		printf("end reached\n");
		isrcstart[it][itdest]=myrange;
	    }else{
		while (start < end-1){
		    auto k= (start + end)/2;
		    //		    printf("%d %d %d %d %d\n", it, itdest, start, end, k );
		    if (bstart[k].getsortkey() < nextval){
			start=k;
		    }else{
			end=k;
		    }
		}
		//		printf("final: %d %d %d %d\n", it, itdest, start, end);
		isrcstart[it][itdest]=end;
	    }
	}
	if (nt == 0){
	    showdt(" range determination");
	    
	    // for(auto it=0; it<nt; it++){
	    // 	auto mystartlocal = it*nwork0;
	    // 	auto  bstart =  key[it];
	    // 	//	printf("isrcstart[%d]:", it);
	    // 	for (auto itdest=0; itdest<nt; itdest++){
	    // 	       	    printf(" %d", isrcstart[it][itdest]);
	    // 	    if (itdest!=0){
	    // 	    	printf(" %g %g", bstart[isrcstart[it][itdest]-1].getsortkey(),
	    // 	    	       bstart[isrcstart[it][itdest]].getsortkey()    );
	    // 	    }
	    
	    // 	}
	    // 	   	printf("\n");
	    // }

	}
#pragma omp barrier	
	auto mydestsize =0;
	for(auto itsrc=0; itsrc < nt; itsrc++){
	    mydestsize += isrcstart[itsrc][it+1]-isrcstart[itsrc][it];
	}
	destsize[it]=mydestsize;
#pragma omp barrier	
	if (it==0){
	    deststart[0]=0;
	    for(auto i=1;i<nt;i++) deststart[i] = deststart[i-1]+destsize[i-1];
	    // for (auto i=0;i<nt;i++){
	    // 	printf("dest start, size[%d]= %d, %d\n", i, deststart[i], destsize[i]);
	    // }
	    for(auto i=0;i<nt; i++) {
		if (maxdestsize < destsize[i]) maxdestsize = destsize[i];
	    }
	}
    }

    KeyValuePair localcopy[nt][maxdestsize];
    T bodylocalcopy[nt][maxdestsize];
    //    showdt(" range determination2");
#pragma omp parallel
    {
	auto it=omp_get_thread_num();
	auto mydestsize = destsize[it];
	int idest=0;
	for(auto itsrc=0; itsrc < nt; itsrc++){
	    auto bstart = key[itsrc];
	    for (auto isrc=isrcstart[itsrc][it];
		 isrc < isrcstart[itsrc][it+1]; isrc++){
		//		printf("to local: %d src= %d dest=%d\n",isrc, idest, isrc);
		localcopy[it][idest] = bstart[isrc];
		idest++;
	    }
	}
	//	printf("sort size=%d %d\n", mydestsize,it);
	sort_bodies(localcopy[it], mydestsize);
	//	printf("sort end %d\n", it);

	if (it==0)showdt(" 2nd sort");
	mydestsize = destsize[it];
	auto offset = deststart[it];
	for(auto idest=0; idest < mydestsize; idest++){
	    bodylocalcopy[it][idest]=b[localcopy[it][idest].key];
	}
	if (it==0)    showdt(" Copyforward");
#pragma omp barrier	
	//	auto mydestsize = destsize[it];
	 offset = deststart[it];
	for(auto idest=0; idest < mydestsize; idest++){
	    //	    printf("it, idest, offst=%d %d %d\n", it, idest, offset);
	    b[offset+idest] = bodylocalcopy[it][idest];
	}
    }
    showdt(" Copyback");
}



int main(int argc, char** argv)
{
    std::mt19937 mt;
    
    // 一様実数分布
    // [-1.0, 1.0)の値の範囲で、等確率に実数を生成する
    std::uniform_real_distribution<> dist1(-1.0, 1.0);
    auto n= atoi(argv[1]);
    fprintf(stderr, "n=%d\n", n);
    auto bodies = new Body[n];
    auto b2 = new Body[n];
    auto b3 = new Body[n];
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
    }
    auto t0=GetWtime();
    sort_bodies(bodies, n);
    auto t1=GetWtime();
    samplesort_bodies(b2, n);
    auto t2a=GetWtime();
    sort_bodies(b2, n);
    auto t2=GetWtime();
    samplesort_bodies(b3, n);
    auto t3=GetWtime();
    printf("Time to sort = %g %g %g\n", t1-t0, t2a-t1, t3-t2);
    bool ok=true;
    for (auto i=0; i<n; i++){
	// printf("%d  %ld %ld, %g %g\n",i,  bodies[i].id, b3[i].id,
	//        bodies[i].pos[0], b3[i].pos[0]);
	if (bodies[i].id != b3[i].id){
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
