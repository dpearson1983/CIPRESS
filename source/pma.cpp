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
#include "../include/cic.h"
#include "../include/pma.h"
#include "../include/pods.h"
#include "../include/particle.h"

#define PI 3.14159265358979
#define G 

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
    
    // Perform initial Fourier transform of density field
    fftw_execute(az_forward);
    
    // Get the individual components of the acceleration ready to transform again
    for (int i = 0; i < N.x; ++i) {
        for (int j = 0; j < N.y; ++j) {
            for (int k = 0; k <= N.z/2; ++k) {
                double k_mag = sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k]);
                
                double rho_real = az[(2*k    ) + 2*(N.z/2 + 1)*(j + N.y*i)];
                double rho_imag = az[(2*k + 1) + 2*(N.z/2 + 1)*(j + N.y*i)];
                
                int index_real = (2*k    ) + 2*(N.z/2 + 1)*(j + N.y*i);
                int index_imag = (2*k + 1) + 2*(N.z/2 + 1)*(j + N.y*i);
                
                ax[index_real] = -4.0*PI*G*rho_real*kx[i]/(k_mag*k_mag);
                ax[index_imag] = -4.0*PI*G*rho_imag*kx[i]/(k_mag*k_mag);
                
                ay[index_real] = -4.0*PI*G*rho_real*ky[i]/(k_mag*k_mag);
                ay[index_imag] = -4.0*PI*G*rho_imag*ky[i]/(k_mag*k_mag);
                
                az[index_real] = -4.0*PI*G*rho_real*kz[i]/(k_mag*k_mag);
                az[index_imag] = -4.0*PI*G*rho_imag*kz[i]/(k_mag*k_mag);
            }
        }
    }
    
    fftw_execute(ax_backward);
    fftw_execute(ay_backward);
    fftw_exectue(az_backward);
}
