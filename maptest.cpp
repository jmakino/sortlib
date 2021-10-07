#include<cmath>
#include<vector>
#include<functional>
#include<algorithm>
#include<exception>
#include<stdexcept>
#include<cassert>
#include<typeinfo>
#include<cstdio>
#include<iostream>
#include<cstring>
#include<map>
#include<random>
#include<omp.h>
#include<stdlib.h>

#include "samplesortlib.hpp"
#include "simplemap.hpp"


int main(int argc, char** argv)
{
    std::mt19937 mt;
    auto nmax= atoi(argv[1]);
    std::uniform_real_distribution<> dist1(nmax);
    auto n= atoi(argv[2]);
    printf("range, n= %d %d\n", nmax, n);
    auto keys = new int[n];
    //    auto m = SimpleMapLib::Map<int64_t, int>(n);
    SimpleMapLib::Map<int64_t, int> m;
    
    std::map<int64_t, int> mstd;
    mstd.clear();
    m.resize(n);
    bool used[nmax];
    for(auto i=0; i< nmax; i++) used[i]=false;
    

    for(auto i=0; i< n; i++){
	int64_t test;
	do{
	    test = dist1(mt);
	}while (used[test]);
	keys[i]=test;
	used[test]=true;
    }
    for(auto i=0; i< n; i++){
	//	m.set(keys[i], i);
	m[keys[i]]= i;
	mstd[keys[i]]= i;
    }
    m.dump();
    m.makemap();
    bool passed=true;
    for(auto i=0; i< n; i++){
	int index = m.at(keys[i]);
	auto p = m.find(keys[i]);
	printf("%d: %d : %d %d\n", keys[i], index, p->second, mstd[keys[i]]);
	if (index != i || p->second!= i || mstd[keys[i]] != i){
	    passed = false;
	    printf("Failed\n");
	}
    }
    if (passed){
	printf("simplemap passed\n");
	return 0;
    }else{
	printf("simplemap failed\n");
	return 1;
    }
}

