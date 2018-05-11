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
#include <fftw3.h>
#include <omp.h>
#include "../include/pma.h"
#include "../include/pods.h"
#include "../include/particle.h"

#define PI 3.14159265358979

std::vector<double> fftFreq(int N, double L) {
    std::vector<double> k;
    k.reserve(N);
    double dk = (2.0*PI)/L;
    for (size_t i = 0; i <= N/2; ++i)
        k.push_back(i*dk);
    for (size_t i = N/2 + 1; i < N; ++i)
        k.push_back((i - N)*dk);
    return k;
}

void getCICInfo(double3 pos, int3 N, double3 L, std::vector<int> &indices, 
                   std::vector<double> &weights) {
    int3 ngp = {int(pos.x/L.x), int(pos.y/L.y), int(pos.z/L.z)};
    double3 del_r = {L.x/N.x, L.y/N.y, L.z/N.z};
    double3 r_ngp = {(ngp.x + 0.5)*del_r.x, (ngp.y + 0.5)*del_r.y, (ngp.z + 0.5)*fel_r.z};
    double3 dr = {pos.x - r_ngp.x, pos.y - r_ngp.y, pos.z - r_ngp.z};
    int3 shift = {dr.x/fabs(dr.x), dr.y/fabs(dr.y), dr.z/fabs(dr.z)};
    
    dr.x = fabs(dr.x);
    dr.y = fabs(dr.y);
    dr.z = fabs(dr.z);
    
    double dV = del_r.x*del_r.y*del_r.z;
    
    indices.push_back(ngp.z + N.z*(ngp.y + N.y*ngp.x));
    indices.push_back(ngp.z + N.z*(ngp.y + N.y*(ngp.x + shift.x)));                         // Shift: x
    indices.push_back(ngp.z + N.z*((ngp.y + shift.y) + N.y*ngp.x));                         // Shift: y
    indices.push_back((ngp.z + shfit.z) + N.z*(ngp.y + N.y*ngp.x));                         // Shift: z
    indices.push_back(ngp.z + N.z*((ngp.y + shift.y) + N.y*(ngp.x + shift.x)));             // Shift: x, y
    indices.push_back((ngp.z + shift.z) + N.z*(ngp.y + N.y*(ngp.x + shift.x)));             // Shift: x, z
    indices.push_back((ngp.z + shift.z) + N.z*((ngp.y + shift.y) + N.y*ngp.x));             // Shift: y, z
    indices.push_back((ngp.z + shift.z) + N.z*((ngp.y + shift.y) + N.y*(ngp.x + shift.x))); // Shift: x, y, z
    
    weights.push_back(((del_r.x - dr.x)*(del_r.y - dr.y)*(del_r.z - dr.z))/dV);
    weights.push_back((dr.x*(del_r.y - dr.y)*(del_r.z - dr.z))/dV);
    weights.push_back(((del_r.x - dr.x)*dr.y*(del_r.z - dr.z))/dV);
    weights.push_back(((del_r.x - dr.x)*(del_r.y - dr.y)*dr.z)/dV);
    weights.push_back((dr.x*dr.y*(del_r.z - dr.z))/dV);
    weights.push_back((dr.x*(del_r.y - dr.y)*dr.z)/dV);
    weights.push_back(((del_r.x - dr.x)*dr.y*dr.z)/dV);
    weights.push_back((dr.x*dr.y*dr.z)/dV);
}

std::vector<double3> particleMeshAcceleration(std::vector<pariticle> &particles, int3 N, double3 L) {
    // Get the frequency vector components
    std::vector<double> kx = fftFreq(N.x, L.x);
    std::vector<double> ky = fftFreq(N.y, L.y);
    std::vector<double> kz = fftFreq(N.z, L.z);
    
    // Set up arrays for the in-place FFTs
    std::vector<double> ax(N.x*N.y*2*(N.z/2 + 1));
    std::vector<double> ay(N.x*N.y*2*(N.z/2 + 1));
    std::vector<double> az(N.x*N.y*2*(N.z/2 + 1));
    
    double dV = (L.x/N.x)*(L.y/N.y)*(L.z/N.z);
    
    // Set up the FFT plans
    fftw_init_threads();
    fftw_import_wisdom_from_filename("fftwWisdom.dat");
    fftw_plan_with_nthreads(omp_get_max_threads());
    fftw_plan ax_forward = fftw_plan_dft_r2c_3d(N.x, N.y, N.z, ax.data(), ax.data(), FFTW_MEASURE);
    fftw_plan ay_forward = fftw_plan_dft_r2c_3d(N.x, N.y, N.z, ay.data(), ay.data(), FFTW_MEASURE);
    fftw_plan az_forward = fftw_plan_dft_r2c_3d(N.x, N.y, N.z, az.data(), az.data(), FFTW_MEASURE);
    fftw_plan ax_backward = fftw_plan_dft_c2r_3d(N.x, N.y, N.z, ax.data(), ax.data(), FFTW_MEASURE);
    fftw_plan ay_backward = fftw_plan_dft_c2r_3d(N.x, N.y, N.z, ay.data(), ay.data(), FFTW_MEASURE);
    fftw_plan az_backward = fftw_plan_dft_c2r_3d(N.x, N.y, N.z, az.data(), az.data(), FFTW_MEASURE);
    fftw_export_wisdom_to_filename("fftwWisdom.dat");
    
    // To keep memory requirements at least somewhat lower, use the az array to calculate the density
    for (std::vector<particle>::iterator it = particles.begin(); it != particles.end; ++it) {
        double3 pos = *it.get_position();
        std::vector<int> indices;
        std::vector<double> weights;
        getCICInfo(pos, N, L, indices, weights);
        
        for (int i = 0; i < 8; ++i)
            az[indices[i]] += weights[i]*(*it.get_mass()/dV);
    }
    
    // TODO: Fourier transform stuff
}
