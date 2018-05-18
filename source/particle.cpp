#include "../particle.h"

double3 particle::get_acceleration(const std::vector<double> &phi, double3 r) {
    

// Implementation of a kick-drift-kick leap-frog integrator.
void particle::update(const std::vector<double> &phi, double dt) {
    double3 acceleration = get_acceleration(phi, this.pos);
    if (!this.firstStep) {
        this.vel.x += this.vel.x + 0.5*acceleration.x*dt;
        this.vel.y += this.vel.y + 0.5*acceleration.y*dt;
        this.vel.z += this.vel.z + 0.5*acceleration.z*dt;
        
        this.pos.x += this.vel.x*dt + 0.5*acceleration.x*dt*dt;
        this.pos.y += this.vel.y*dt + 0.5*acceleration.y*dt*dt;
        this.pos.z += this.vel.z*dt + 0.5*acceleration.z*dt*dt;
        
        this.vel.x += this.vel.x + 0.5*acceleration.x*dt;
        this.vel.y += this.vel.y + 0.5*acceleration.y*dt;
        this.vel.z += this.vel.z + 0.5*acceleration.z*dt;
    } else {
        this.pos.x += this.vel.x*dt + 0.5*acceleration.x*dt*dt;
        this.pos.y += this.vel.y*dt + 0.5*acceleration.y*dt*dt;
        this.pos.z += this.vel.z*dt + 0.5*acceleration.z*dt*dt;
        
        this.vel.x += this.vel.x + 0.5*acceleration.x*dt;
        this.vel.y += this.vel.y + 0.5*acceleration.y*dt;
        this.vel.z += this.vel.z + 0.5*acceleration.z*dt;
        
        this.firstStep = false;
    }
}

double3 particle::get_position() {
    return this.pos;
}

double3 particle::get_velocity() {
    return this.vel;
}

double particle::get_mass() {
    return this.mass;
}

size_t particle::get_ID() {
    return this.ID;
}
