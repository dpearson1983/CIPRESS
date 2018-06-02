/* pma.cpp
 * David W. Pearson
 * 11 May 2018
 * 
 * This code will calculate the gravitational field of a large collection of particles through
 * the use of Fourier transforms. To start, the particles masses are binned on a grid. Then the 
 * value at each grid cell will be divided by the volume of that cell to yield a density field.
 * The density field is then Fourier transformed, multiplied by the frequency vector along with
 * the constants (-4*pi*G), divided by the magnitude of the frequency squared, and then inverse
 * transformed (this is actually three different inverse transforms, one for each vector 
 * component). The functions return a vector of double3 type storing a.x, a.y, and a.z on a
 * grid.
 * 
 * Also contained here is the code to interpolate the acceleration to a specific point based on
 * the values of the set of nearest grid points.
 */


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <fftw3.h>
#include <omp.h>
#include "../include/cic.h"
#include "../include/pma.h"
#include "../include/pods.h"
#include "../include/particle.h"

#ifndef PI
#define PI 3.14159265358979
#endif

#ifndef G
#define G 43215.3390466 // In Gadget2 style units.
#endif

std::vector<double> fftFreq(const int N, const double L) {
    std::vector<double> k;
    k.reserve(N);
    double dk = (2.0*PI)/L;
    for (size_t i = 0; i <= N/2; ++i)
        k.push_back(i*dk);
    for (size_t i = N/2 + 1; i < N; ++i)
        k.push_back((i - N)*dk);
    return k;
}

void particleMeshAcceleration(std::vector<double> &rho, std::vector<fftw_complex> &phi, 
                              const std::vector<particle> &particles, const std::vector<double> &kx, 
                              const std::vector<double> &ky, const std::vector<double> &kz, const int3 &N,
                              const double3 &L, std::string wisdomFile) {
    double dV = (L.x/N.x)*(L.y/N.y)*(L.z/N.z);
    size_t N_tot = N.x*N.y*N.z;
    
    // Set up the FFT plans
    fftw_init_threads();
    fftw_import_wisdom_from_filename(wisdomFile.c_str());
    fftw_plan_with_nthreads(omp_get_max_threads());
    fftw_plan forward = fftw_plan_dft_r2c_3d(N.x, N.y, N.z, rho.data(), phi.data(), FFTW_MEASURE);
    fftw_plan backward = fftw_plan_dft_c2r_3d(N.x, N.y, N.z, phi.data(), rho.data(), FFTW_MEASURE);
    fftw_export_wisdom_to_filename(wisdomFile.c_str());
    
    // Bin the particles using cloud-in-cell to get density field
    for (auto it = particles.begin(); it != particles.end(); ++it) {
        particle p = *it;
        double3 pos = p.get_position();
        std::vector<size_t> indices;
        std::vector<double> weights;
        getCICInfo(pos, N, L, indices, weights);
        
        for (int i = 0; i < 8; ++i)
            rho[indices[i]] += weights[i]*(p.get_mass()/dV);
    }
    
    // Perform initial Fourier transform of density field
    fftw_execute(forward);
    
    // Get the individual components of the acceleration ready to transform again
    for (int i = 0; i < N.x; ++i) {
        for (int j = 0; j < N.y; ++j) {
            for (int k = 0; k <= N.z/2; ++k) {
                double k_mag = sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k]);
                int index = k + (N.z/2 + 1)*(j + N.y*i);
                
                if (k_mag > 0) {                
                    phi[index][0] = -4.0*PI*G/(k_mag*k_mag);
                    phi[index][1] = -4.0*PI*G/(k_mag*k_mag);
                } else {
                    phi[index][0] = 0.0;
                    phi[index][1] = 0.0;
                }
            }
        }
    }
    
    // Perform backward transform to get the potential field
    fftw_execute(backward);
    
    // Normalize
    for (size_t i = 0; i < N_tot; ++i)
        rho[i] /= N_tot;
    
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_cleanup();
}

void particleMeshAcceleration(std::vector<double> &rho, std::vector<fftw_complex> &phi, 
                              particle *particles, const int N_p, const std::vector<double> &kx, 
                              const std::vector<double> &ky, const std::vector<double> &kz, const int3 &N,
                              const double3 &L, std::string wisdomFile) {
    double dV = (L.x/N.x)*(L.y/N.y)*(L.z/N.z);
    size_t N_tot = N.x*N.y*N.z;
    
    // Set up the FFT plans
    fftw_init_threads();
    fftw_import_wisdom_from_filename(wisdomFile.c_str());
    fftw_plan_with_nthreads(omp_get_max_threads());
    fftw_plan forward = fftw_plan_dft_r2c_3d(N.x, N.y, N.z, rho.data(), phi.data(), FFTW_MEASURE);
    fftw_plan backward = fftw_plan_dft_c2r_3d(N.x, N.y, N.z, phi.data(), rho.data(), FFTW_MEASURE);
    fftw_export_wisdom_to_filename(wisdomFile.c_str());
    
    // Bin the particles using cloud-in-cell to get density field
    for (size_t i = 0; i < N_p; ++i) {
        double3 pos = particles[i].get_position();
        std::vector<size_t> indices;
        std::vector<double> weights;
        getCICInfo(pos, N, L, indices, weights);
        
        for (int i = 0; i < 8; ++i)
            rho[indices[i]] += (weights[i]*particles[i].get_mass())/dV;
    }
    
    // Perform initial Fourier transform of density field
    fftw_execute(forward);
    
    // Get the individual components of the acceleration ready to transform again
    for (int i = 0; i < N.x; ++i) {
        for (int j = 0; j < N.y; ++j) {
            for (int k = 0; k <= N.z/2; ++k) {
                double k_mag = sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k]);
                int index = k + (N.z/2 + 1)*(j + N.y*i);
                
                if (k_mag > 0) {                
                    phi[index][0] = -4.0*PI*G/(k_mag*k_mag);
                    phi[index][1] = -4.0*PI*G/(k_mag*k_mag);
                } else {
                    phi[index][0] = 0.0;
                    phi[index][1] = 0.0;
                }
            }
        }
    }
    
    // Perform backward transform to get the potential field
    fftw_execute(backward);
    
    // Normalize
    for (size_t i = 0; i < N_tot; ++i)
        rho[i] /= N_tot;
    
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_cleanup();
}

// This function imposes periodic boundary conditions for an arbitrary shift
size_t periodic(size_t i, int shift, int N) {
    size_t ip = i + shift % N;
    if (ip > 2*N) {
        ip = N - (std::numeric_limits<size_t>::max() - ip) - 1;
    } else if (ip >= N) {
        ip = ip % N;
    }
    return ip;
}

double3 gradient(size_t4 index, const std::vector<double> &phi, const int3 N, const double3 dr) {
    double3 a = {0.0, 0.0, 0.0};
    
    a.x = (phi[index.z + N.z*(index.y + N.y*periodic(index.x, -2, N.x))] -
           8.0*phi[index.z + N.z*(index.y + N.y*periodic(index.x, -1, N.x))] +
           8.0*phi[index.z + N.z*(index.y + N.y*periodic(index.x, 1, N.x))] -
           phi[index.z + N.z*(index.y + N.y*periodic(index.x, 2, N.x))])/(12.0*dr.x);
    a.y = (phi[index.z + N.z*(periodic(index.y, -2, N.y) + N.y*index.x)] -
           8.0*phi[index.z + N.z*(periodic(index.y, -1, N.y) + N.y*index.x)] +
           8.0*phi[index.z + N.z*(periodic(index.y, 1, N.y) + N.y*index.x)] -
           phi[index.z + N.z*(periodic(index.y, 2, N.y) + N.y*index.x)])/(12.0*dr.y);
    a.z = (phi[periodic(index.z, -2, N.z) + N.z*(index.y + N.y*index.x)] -
           8.0*phi[periodic(index.z, -1, N.z) + N.z*(index.y + N.y*index.x)] +
           8.0*phi[periodic(index.z, 1, N.z) + N.z*(index.y + N.y*index.x)] -
           phi[periodic(index.z, 2, N.z) + N.z*(index.y + N.y*index.x)])/(12.0*dr.x);
           
    return a;
}

// Calculates the acceleration of the particle at position by taking the gradient of the gravitational
// potential using the 5 point finite difference method.
double3 getParticleAcceleration(double3 position, const std::vector<double> &phi, const int3 &N, 
                                const double3 &L) {
    std::vector<size_t4> indices;
    std::vector<double> weights;
    getCICInfo(position, N, L, indices, weights);
    double3 a = {0.0, 0.0, 0.0};
    
    double3 dr = {L.x/N.x, L.y/N.y, L.z/N.z};
    for (int i = 0; i < 8; ++i) {
        double3 a_grid = gradient(indices[i], phi, N, dr);
        a.x += weights[i]*a_grid.x;
        a.y += weights[i]*a_grid.y;
        a.z += weights[i]*a_grid.z;
    }
    return a;
}

// NOTE: This function is for testing purposes only. It will be very slow but very precise. It is
//       intended to be used to compare the results of the PM methods above with an "exact" result using
//       the particle-particle method.
double3 getParticleAcceleration(double3 position, const std::vector<particle> &particles) {
    double3 a = {0.0, 0.0, 0.0};
    for (auto it = particles.begin(); it != particles.end(); ++it) {
        particle p = *it;
        double3 r = p.get_position() - position;
        double r_mag = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
        
        a.x += G*(p.get_mass())*r.x/(r_mag*r_mag*r_mag);
        a.y += G*(p.get_mass())*r.y/(r_mag*r_mag*r_mag);
        a.z += G*(p.get_mass())*r.z/(r_mag*r_mag*r_mag);
    }
    return a;
}
