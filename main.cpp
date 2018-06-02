#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <fftw3.h>
#include <omp.h>
#include <harppi.h>
#include <include/cic.h>
#include <include/pma.h>
#include <include/pods.h>
#include <include/particle.h>
#include <include/gadgetReader.h>

std::vector<particle> getParticles(std::string particleFile, double3 &L) {
    std::vector<std::vector<float3>> pos(6);
    std::vector<std::vector<float3>> vel(6);
    std::vector<std::vector<int>> ids(6);
    gadget_header header;
    
    read_gadget_snapshot(particleFile, header, pos, vel, ids);
    
    print_gadget_header(header);
    
    L.x = L.y = L.z = header.boxSize;
    
    std::vector<particle> particles;
    particles.reserve(header.N_p[1]);
    for (int i = 0; i < header.N_p[1]; ++i) {
        particle p(ids[1][i], pos[1][i], vel[1][i], header.M_p[1]);
        particles.push_back(p);
    }
    
    return particles;
}

int main(int argc, char *argv[]) {
    parameters p(argv[1]);
    
    int3 N = {p.geti("Nx"), p.geti("Ny"), p.geti("Nz")};
    double3 L;
    
    std::vector<double> rho(N.x*N.y*N.z);
    std::vector<fftw_complex> phi(N.x*N.y*(N.z/2 + 1));
    
    std::vector<particle> particles = getParticles(p.gets("snapshotFile"), L);
    
    std::vector<double> kx = fftFreq(N.x, L.x);
    std::vector<double> ky = fftFreq(N.y, L.y);
    std::vector<double> kz = fftFreq(N.z, L.z);
    
    std::random_device seeder;
    std::mt19937_64 gen(seeder());
    std::uniform_real_distribution<double> dist (0.0, 1.0);
    
    particleMeshAcceleration(rho, phi, particles, kx, ky, kz, N, L, p.gets("wisdomFile"));
    
    std::ofstream fout(p.gets("outFile"));
    fout.precision(15);
    for (int i = 0; i < p.geti("numTests"); ++i) {
        double3 pos = {dist(gen)*L.x, dist(gen)*L.y, dist(gen)*L.z};
        
        double3 a_pp = getParticleAcceleration(pos, particles);
        double3 a_phi = getParticleAcceleration(pos, rho, N, L);
        
        double3 diff = a_pp - a_phi;
        
        fout << a_pp.x << " " << a_pp.y << " " << a_pp.z << " " << a_phi.x << " " << a_phi.y << " ";
        fout << a_phi.z << " " << diff.x << " " << diff.y << " " << diff.z << "\n";
    }
    fout.close();
    
    return 0;
}
