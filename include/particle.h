#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "pods.h"

class particle{
    double3 pos, vel;
    double mass;
    size_t ID;
    
    public:
        void update(double3 acceleration, double dt);
        
        double3 get_position();
        
        double3 get_velocity();
        
        double get_mass();
        
        size_t get_ID();
};

#endif
