#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <vector>
#include "pods.h"

// Class for storing relevant information for particles in N-body simulations. Each particle will need
// roughly 64 bytes of memory for storage.
//      128^3 particles =    0.125 GiB =    0.134217728 GB
//      256^3 particles =    1     GiB =    1.073741824 GB
//      512^3 particles =    8     GiB =    8.589934592 GB
//     1024^3 particles =   64     GiB =   68.719476736 GB
//     2048^3 particles =  512     GiB =  549.755813888 GB
//     4096^3 particles = 4096     GiB = 4398.0465111   GB
class particle{
    double3 pos, vel; // 48 bytes
    double mass; // 8 bytes
    size_t ID; // 8 bytes
    
    public:
        particle(const int &id, const double3 &position, const double3 &velocity, const double &m);
        
        particle(const int &id, const float3 &posistion, const float3 &velocity, const double &m);
        
        void firstStep(const std::vector<double> &phi, const double dt, const int3 &N, const double3 &L);
        
        void update(const std::vector<double> &phi, const double dt, const int3 &N, const double3 &L);
        
        double3 get_position();
        
        double3 get_velocity();
        
        double get_mass();
        
        size_t get_ID();
};

#endif
