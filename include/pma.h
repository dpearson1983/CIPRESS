#ifndef _PMA_H_
#define _PMA_H_

#include <vector>
#include <string>
#include <fftw3.h>
#include "pods.h"
#include "particle.h"

std::vector<double> fftFreq(const int N, const double L);

void particleMeshAcceleration(std::vector<double> &rho, std::vector<fftw_complex> &phi, 
                              const std::vector<particle> &particles, const std::vector<double> &kx, 
                              const std::vector<double> &ky, const std::vector<double> &kz, const int3 &N,
                              const double3 &L, std::string wisdomFile = "fftwWisdom.dat", 
                              bool verb = false);

void particleMeshAcceleration(std::vector<double> &rho, std::vector<fftw_complex> &phi, 
                              particle *particles, const int N_p, const std::vector<double> &kx, 
                              const std::vector<double> &ky, const std::vector<double> &kz, const int3 &N,
                              const double3 &L, std::string wisdomFile = "fftwWisdom.dat", 
                              bool verb = false);

double3 getParticleAcceleration(double3 position, const std::vector<double> &phi, const int3 &N, 
                                const double3 &L);

double3 getParticleAcceleration(double3 position, const std::vector<particle> &particles);

#endif
