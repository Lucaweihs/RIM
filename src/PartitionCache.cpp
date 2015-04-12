#ifndef RIM_PartitionCache
#define RIM_PartitionCache

#include<stdio.h>
#include<stdlib.h>
#include<cstdlib>
#include<algorithm>
#include<math.h>
#include"List.cpp"
//#include<boost/unordered_map.hpp>
//#typedef boost::unordered_map<>

namespace RIM {
  class PartitionCache {
    private:
      int* partMat;
      int maxN;
      int maxLimit;
      int maxParts;
      
      bool validInput(int n, int limit, int parts) {
        if(n > maxN || limit > maxLimit || parts > maxParts ||
           n < 0 || limit < 0 || parts < 0) {
             return(false);
        }
        return(true);
      }
      void setIndex(int val, int n, int limit, int parts) {
        if(!validInput(n, limit, parts)) {
          std::exit(1);
        }
        partMat[n*(maxLimit+1)*(maxLimit+1) + limit*(maxParts+1) + parts] = val;
      }
      void addToIndex(int toAdd, int n, int limit, int parts) {
        if(!validInput(n, limit, parts)) {
          std::exit(1);
        }
        partMat[n*(maxLimit+1)*(maxLimit+1) + limit*(maxParts+1) + parts] += toAdd;
      }
      
    public:
      PartitionCache(int n, int limit, int parts) {
        maxN = n;
        maxParts = parts;
        maxLimit = limit;
        partMat = (int*) malloc(((n+1)*(parts+1)*(limit+1)+1)*sizeof(int));
        for(int i=0; i<(n+1)*(parts+1)*(limit+1)+1; i++) {
            partMat[i] = -1;
        }
      }
      
      int atIndex(int n, int limit, int parts) {
        if(!validInput(n, limit, parts)) {
          std::exit(1);
        }
        return(partMat[n*(maxLimit+1)*(maxLimit+1) + limit*(maxParts+1) + parts]);
      }
      
      bool hasBeenComputed(int n, int limit, int parts) {
        if(!validInput(n, limit, parts)) {
          std::exit(1);
        }
        return(atIndex(n, limit, parts) != -1);
      }
      
      int countPartitions(int n, int limit, int parts) {
        if(!validInput(n, limit, parts)) {
          std::exit(1);
        }
        
        if(hasBeenComputed(n,limit,parts)) {
          return(atIndex(n,limit,parts));
        } else if(n == 0) {
          setIndex(1,n,limit,parts);
          return(1);
        } else if(n > limit*parts) {
          setIndex(0,n,limit,parts);
          return(0);
        }  else {
          setIndex(0,n,limit,parts);
          for(int newLimit=1; newLimit <= std::min(limit, n); newLimit++) {
            addToIndex(this->countPartitions(n-newLimit, newLimit, parts-1), n, limit, parts);
          }
          return(atIndex(n, limit, parts));
        }
      }
      
      RIM::List<int>* randomPartition(double rand, int n, int limit, int parts) {
        RIM::List<int>* l = new RIM::List<int>();
        int total = countPartitions(n, limit, parts);
        if(total == 0) {
          return(l);
        } else if(rand <= 0 || rand >= 1) {
          std::exit(1);
        }
        
        int which = ceil(rand*total);
        int i, count;
        
        while(n != 0) {
          for(i=1; i<=std::min(limit, n); i++) {
            count = countPartitions(n-i, i, parts-1);
            if(count >= which) {
              break;
            }
            which -= count;
          }
          l->appendValue(i);
          n -= i;
          limit = i;
          parts -= 1;
        }
        return(l);
      }
  };
}

#endif

/*
List* randomPartition(int n, int limit, int maxParts) {
  
}

int countPartitions(int* countMat, int n, int limit, int maxParts, int currentN, int currentLimit, int currentMaxParts) {
  int ind = currentN*(limit*maxParts) + currentLimit*maxParts + currentMaxParts;
  if(countMat[ind] != 0) {
    return(countMat[ind]);
  }
  else if(currentN == 0) {
    countMat[ind] = 1;
    return(1);
  }
  else if(currentN > currentLimit*currentMaxParts) {
    return(0);
  } 
  else {
    for(int newLimit=1; newLimit <= min(currentLimit, currentN); newLimit++) {
      countMat[ind] += countPartitions(countMat, n, limit, maxParts, currentN-newLimit, newLimit, currentMaxParts-1);
    }
    return(countMat[ind]);
  }
}

*/