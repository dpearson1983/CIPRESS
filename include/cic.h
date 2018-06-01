#ifndef _CIC_H_
#define _CIC_H_

#include "pods.h"

void getCICInfo(double3 pos, const int3 &N, const double3 &L, std::vector<size_t> &indices, 
                std::vector<double> &weights);

void getCICInfo(double3 pos, const int3 &N, const double3 &L, std::vector<size_t4> &indices, 
                std::vector<double> &weights);

#endif
