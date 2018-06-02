#include <vector>
#include "../include/particle.h"
#include "../include/cic.h"
#include "../include/pma.h"
#include "../include/pods.h"

particle::particle(const int &id, const double3 &position, const double3 &velocity, const double &m) {
    particle::pos.x = position.x;
    particle::pos.y = position.y;
    particle::pos.z = position.z;
    particle::vel.x = velocity.x;
    particle::vel.y = velocity.y;
    particle::vel.z = velocity.z;
    particle::mass = m;
    particle::ID = id;
}

particle::particle(const int &id, const float3 &position, const float3 &velocity, const double &m) {
    particle::pos.x = position.x;
    particle::pos.y = position.y;
    particle::pos.z = position.z;
    particle::vel.x = velocity.x;
    particle::vel.y = velocity.y;
    particle::vel.z = velocity.z;
    particle::mass = m;
    particle::ID = id;
}

// Do a drift-kick for the first time step.
void particle::firstStep(const std::vector<double> &phi, const double dt, const int3 &N, const double3 &L) {
    double3 acceleration = getParticleAcceleration(particle::pos, phi, N, L);
    
    particle::pos.x += particle::vel.x*dt + 0.5*acceleration.x*dt*dt;
    particle::pos.y += particle::vel.y*dt + 0.5*acceleration.y*dt*dt;
    particle::pos.z += particle::vel.z*dt + 0.5*acceleration.z*dt*dt;
    
    particle::vel.x += particle::vel.x + 0.5*acceleration.x*dt;
    particle::vel.y += particle::vel.y + 0.5*acceleration.y*dt;
    particle::vel.z += particle::vel.z + 0.5*acceleration.z*dt;
}

// Implementation of a kick-drift-kick leap-frog integrator.
void particle::update(const std::vector<double> &phi, const double dt, const int3 &N, const double3 &L) {
    double3 acceleration = getParticleAcceleration(particle::pos, phi, N, L);
    
    particle::vel.x += particle::vel.x + 0.5*acceleration.x*dt;
    particle::vel.y += particle::vel.y + 0.5*acceleration.y*dt;
    particle::vel.z += particle::vel.z + 0.5*acceleration.z*dt;
    
    particle::pos.x += particle::vel.x*dt + 0.5*acceleration.x*dt*dt;
    particle::pos.y += particle::vel.y*dt + 0.5*acceleration.y*dt*dt;
    particle::pos.z += particle::vel.z*dt + 0.5*acceleration.z*dt*dt;
    
    particle::vel.x += particle::vel.x + 0.5*acceleration.x*dt;
    particle::vel.y += particle::vel.y + 0.5*acceleration.y*dt;
    particle::vel.z += particle::vel.z + 0.5*acceleration.z*dt;
}

double3 particle::get_position() {
    return particle::pos;
}

double3 particle::get_velocity() {
    return particle::vel;
}

double particle::get_mass() {
    return particle::mass;
}

size_t particle::get_ID() {
    return particle::ID;
}
