#include <vector>
#include <cmath>
#include "../include/pods.h"

// TODO: 
//  -- Add check for when particle is at center of grid cell?
void getCICInfo(const double3 &pos, const int3 &N, const double3 &L, std::vector<size_t> &indices, 
                   std::vector<double> &weights) {
    double3 del_r = {L.x/N.x, L.y/N.y, L.z/N.z};
    int3 ngp = {int(pos.x/del_r.x), int(pos.y/del_r.y), int(pos.z/del_r.z)};
    double3 r_ngp = {(ngp.x + 0.5)*del_r.x, (ngp.y + 0.5)*del_r.y, (ngp.z + 0.5)*del_r.z};
    double3 dr = {pos.x - r_ngp.x, pos.y - r_ngp.y, pos.z - r_ngp.z};
    int3 shift = {int(dr.x/fabs(dr.x)), int(dr.y/fabs(dr.y)), int(dr.z/fabs(dr.z))};
    
    dr.x = fabs(dr.x);
    dr.y = fabs(dr.y);
    dr.z = fabs(dr.z);
    
    double dV = del_r.x*del_r.y*del_r.z;
    
    if (ngp.x + shift.x == -1) shift.x = N.x - 1;
    if (ngp.y + shift.y == -1) shift.y = N.y - 1;
    if (ngp.z + shift.z == -1) shift.z = N.z - 1;
    
    if (ngp.x + shift.x == N.x) shift.x = 1 - N.x;
    if (ngp.y + shift.y == N.y) shift.y = 1 - N.y;
    if (ngp.z + shift.z == N.z) shift.z = 1 - N.z;
    
    indices.push_back(size_t(ngp.z + N.z*(ngp.y + N.y*ngp.x)));
    indices.push_back(size_t(ngp.z + N.z*(ngp.y + N.y*(ngp.x + shift.x))));                 // Shift: x
    indices.push_back(size_t(ngp.z + N.z*((ngp.y + shift.y) + N.y*ngp.x)));                 // Shift: y
    indices.push_back(size_t((ngp.z + shift.z) + N.z*(ngp.y + N.y*ngp.x)));                 // Shift: z
    indices.push_back(size_t(ngp.z + N.z*((ngp.y + shift.y) + N.y*(ngp.x + shift.x))));     // Shift: x, y
    indices.push_back(size_t((ngp.z + shift.z) + N.z*(ngp.y + N.y*(ngp.x + shift.x))));     // Shift: x, z
    indices.push_back(size_t((ngp.z + shift.z) + N.z*((ngp.y + shift.y) + N.y*ngp.x)));     // Shift: y, z
    indices.push_back(size_t((ngp.z + shift.z) + N.z*((ngp.y + shift.y) + N.y*(ngp.x + shift.x)))); // Shift: x, y, z
    
    weights.push_back(((del_r.x - dr.x)*(del_r.y - dr.y)*(del_r.z - dr.z))/dV);
    weights.push_back((dr.x*(del_r.y - dr.y)*(del_r.z - dr.z))/dV);
    weights.push_back(((del_r.x - dr.x)*dr.y*(del_r.z - dr.z))/dV);
    weights.push_back(((del_r.x - dr.x)*(del_r.y - dr.y)*dr.z)/dV);
    weights.push_back((dr.x*dr.y*(del_r.z - dr.z))/dV);
    weights.push_back((dr.x*(del_r.y - dr.y)*dr.z)/dV);
    weights.push_back(((del_r.x - dr.x)*dr.y*dr.z)/dV);
    weights.push_back((dr.x*dr.y*dr.z)/dV);
}

void getCICInfo(double3 pos, const int3 &N, const double3 &L, std::vector<size_t4> &indices, 
                   std::vector<double> &weights) {
    double3 del_r = {L.x/N.x, L.y/N.y, L.z/N.z};
    int3 ngp = {int(pos.x/del_r.x), int(pos.y/del_r.y), int(pos.z/del_r.z)};
    double3 r_ngp = {(ngp.x + 0.5)*del_r.x, (ngp.y + 0.5)*del_r.y, (ngp.z + 0.5)*del_r.z};
    double3 dr = {pos.x - r_ngp.x, pos.y - r_ngp.y, pos.z - r_ngp.z};
    int3 shift = {int(dr.x/fabs(dr.x)), int(dr.y/fabs(dr.y)), int(dr.z/fabs(dr.z))};
    
    dr.x = fabs(dr.x);
    dr.y = fabs(dr.y);
    dr.z = fabs(dr.z);
    
    double dV = del_r.x*del_r.y*del_r.z;
    
    if (ngp.x + shift.x == -1) shift.x = N.x - 1;
    if (ngp.y + shift.y == -1) shift.y = N.y - 1;
    if (ngp.z + shift.z == -1) shift.z = N.z - 1;
    
    if (ngp.x + shift.x == N.x) shift.x = 1 - N.x;
    if (ngp.y + shift.y == N.y) shift.y = 1 - N.y;
    if (ngp.z + shift.z == N.z) shift.z = 1 - N.z;
    
    size_t4 index1 = {size_t(ngp.x), size_t(ngp.y), size_t(ngp.z), size_t(ngp.z + N.z*(ngp.y + N.y*ngp.x))};
    size_t4 index2 = {size_t(ngp.x + shift.x), size_t(ngp.y), size_t(ngp.z), size_t(ngp.z + N.z*(ngp.y + N.y*(ngp.x + shift.x)))};
    size_t4 index3 = {size_t(ngp.x), size_t(ngp.y + shift.y), size_t(ngp.z), size_t(ngp.z + N.z*((ngp.y + shift.y) + N.y*ngp.x))};
    size_t4 index4 = {size_t(ngp.x), size_t(ngp.y), size_t(ngp.z + shift.z), size_t((ngp.z + shift.z) + N.z*(ngp.y + N.y*ngp.x))};
    size_t4 index5 = {size_t(ngp.x + shift.x), size_t(ngp.y + shift.y), size_t(ngp.z), size_t(ngp.z + N.z*((ngp.y + shift.y) + N.y*(ngp.x + shift.x)))};
    size_t4 index6 = {size_t(ngp.x + shift.x), size_t(ngp.y), size_t(ngp.z + shift.z), size_t((ngp.z + shift.z) + N.z*(ngp.y + N.y*(ngp.x + shift.x)))};
    size_t4 index7 = {size_t(ngp.x), size_t(ngp.y + shift.y), size_t(ngp.z + shift.z), size_t((ngp.z + shift.z) + N.z*((ngp.y + shift.y) + N.y*ngp.x))};
    size_t4 index8 = {size_t(ngp.x + shift.x), size_t(ngp.y + shift.y), size_t(ngp.z + shift.z), size_t((ngp.z + shift.z) + N.z*((ngp.y + shift.y) + N.y*(ngp.x + shift.x)))};
    
    indices.push_back(index1);
    indices.push_back(index2);                         // Shift: x
    indices.push_back(index3);                         // Shift: y
    indices.push_back(index4);                         // Shift: z
    indices.push_back(index5);             // Shift: x, y
    indices.push_back(index6);             // Shift: x, z
    indices.push_back(index7);             // Shift: y, z
    indices.push_back(index8); // Shift: x, y, z
    
    weights.push_back(((del_r.x - dr.x)*(del_r.y - dr.y)*(del_r.z - dr.z))/dV);
    weights.push_back((dr.x*(del_r.y - dr.y)*(del_r.z - dr.z))/dV);
    weights.push_back(((del_r.x - dr.x)*dr.y*(del_r.z - dr.z))/dV);
    weights.push_back(((del_r.x - dr.x)*(del_r.y - dr.y)*dr.z)/dV);
    weights.push_back((dr.x*dr.y*(del_r.z - dr.z))/dV);
    weights.push_back((dr.x*(del_r.y - dr.y)*dr.z)/dV);
    weights.push_back(((del_r.x - dr.x)*dr.y*dr.z)/dV);
    weights.push_back((dr.x*dr.y*dr.z)/dV);
}
