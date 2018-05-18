#ifndef _PMA_H_
#define _PMA_H_

#include <vector>
#include <string>
#include "pods.h"
#include "particle.h"

void particleMeshAcceleration(std::vector<double> &phi, std::vector<particle> &particles, int3 N, double3 L);

void particleMeshAcceleration(std::vector<double> &phi, std::string particle_file, int3 N, double3 L);

double3 getParticleAcceleration(double3 position);

#endif
