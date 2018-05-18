#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "pods.h"

// Class for storing relevant information for particles in N-body simulations. Each particle will need
// roughly 68 bytes of memory for storage.
//      128^3 particles =    0.1328125 GiB =    0.142606336 GB
//      256^3 particles =    1.0625    GiB =    1.14085069  GB
//      512^3 particles =    8.5       GiB =    9.12681     GB
//     1024^3 particles =   68         GiB =   73.0144      GB
//     2048^3 particles =  544         GiB =  584.116       GB
//     4096^3 particles = 4352         GiB = 4672.924       GB
class particle{
    double3 pos, vel; // 48 bytes
    double mass; // 8 bytes
    size_t ID; // 8 bytes
    bool firstStep = true; // 4 bytes
    
    double3 get_acceleration(const std::vector<double> &phi, double3 r);
    
    public:
        void update(const std::vector<double> &phi, double dt);
        
        double3 get_position();
        
        double3 get_velocity();
        
        double get_mass();
        
        size_t get_ID();
};

#endif
