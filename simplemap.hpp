#pragma once
#include <unistd.h>
//#include "samplesortlib.hpp"
//
// simplemap.hpp
//
// a simple (but OMP-friendly) replacement of std::map
//
// provide:
//
//#define SIMPLEMAP_RANGE_CHECK	    
namespace SimpleMapLib{

    const int slack=1024;
    template <typename KeyT, typename ValT>
    class Map{
    public:

	class KeyValuePair{
	public:
	    KeyT first;
	    ValT second;
	};

	class KeyAddressPair{
	public:
	    KeyT first;
	    KeyValuePair * addr;
	    ValT operator =(ValT v){
		addr[v].first=first;
		addr[v].second=v;
		return v;
	    }
	    KeyAddressPair(KeyT k, KeyValuePair* a){
		first=k;
		addr = a;
	    }
	    KeyAddressPair(){}
	};
	
	int32_t size;
	int32_t memsize;
	KeyValuePair * p;
	KeyAddressPair *passignbuf;
	bool sorted;
	void initialize(ValT s){
	    if (s>0){
		p = new KeyValuePair[s+slack];
		size=s;
		memsize = s+slack;
#ifdef _OPENMP
		passignbuf = new KeyAddressPair[omp_get_max_threads()];
#else
		passignbuf = new KeyAddressPair[1];
#endif		
		
	    }else{
		p= NULL;
		size=0;
		memsize=0;
	    }
	    sorted = false;
	}
    
	Map(ValT s){ initialize(s);	}
	Map(){ initialize(0);	}
	~Map(){
	    if (p!=NULL)  delete[] p;
	}
	void resize(ValT s){
	    sorted = false;
	    size=s;
	    if (s >memsize){
		if (p!=NULL) delete[] p;
		p = new KeyValuePair[s+slack];
		memsize = s+slack;
	    }
	}
	void clear(){
	    sorted=false;
	    size=0;
	}
		
		
	void makemap()
	{
#if 0	    
	    printf("data before sort\n");
	    for(auto i=0;i<size; i++){
		printf("it=%d i= %d  key=%d val=%d\n",
		       omp_get_thread_num(),
		       i, (int) p[i].first, (int) p[i].second);
	    }
#endif	    
		
	    SampleSortLib::samplesort(p, size,[](KeyValuePair & l)
				->auto{return l.first;});
	    sorted=true;
	}
	void set(KeyT  k, ValT v)
	{
#ifdef SIMPLEMAP_RANGE_CHECKE	    
	    if (v >= size){
		printf("simplemap::set failed, too large index:%d size:%d %d\n",
		       v, size, memsize);
		exit(-1);
	    }
#endif	    
	    p[v].first=k;
	    p[v].second=v;
	}

	KeyValuePair* find(KeyT k){
	    if (!sorted)makemap();
		
	    if (k < p[0].first || k > p[size-1].first ){
		return p+size;
	    }
	    ValT low=0;
	    ValT high=size-1;
	    auto mid = (low + high)/2;
	    while (low < high-1){
		if (p[mid].first < k){
		    low=mid;
		}else{
		    high=mid;
		}
		mid = (low + high)/2;
		if (p[mid].first == k){
		    return  p+mid;
		}
	    }
	    if (p[high].first == k){
		return  p+high;
	    }else if (p[low].first == k){
		return  p+low;
	    }
	    return p+(size);

	}
	ValT value(KeyT k){return find(k)->second;}
	ValT at(KeyT k){return find(k)->second;}
	KeyValuePair* end(){return p+size;}

	KeyAddressPair  operator [](KeyT  k) {
	    return KeyAddressPair(k,p);
	}
	//	const ValT  operator [](KeyT  k) {return at(k);}
	void dump()
	{
	    printf("Size: %d\n", size);
	    for(auto i=0;i<size; i++){
		printf("%ld: %d\n", (int64_t)p[i].first, p[i].second);
	    }
	}
    };
}
	    

    
