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
    
    L.x = L.y = L.z = header.boxSize/1000.0;
    
    std::vector<particle> particles;
    particles.reserve(header.N_p[1]);
    for (int i = 0; i < header.N_p[1]; ++i) {
        float3 position = {pos[1][i].x/1000.0, pos[1][i].y/1000.0, pos[1][i].z/1000.0};
        particle p(ids[1][i], position, vel[1][i], header.M_p[1]);
        particles.push_back(p);
    }
    
    return particles;
}

int main(int argc, char *argv[]) {
    parameters p(argv[1]);
    p.print();
    
    int3 N = {p.geti("Nx"), p.geti("Ny"), p.geti("Nz")};
    double3 L;
    
    std::vector<double> rho(N.x*N.y*N.z);
    std::vector<fftw_complex> phi(N.x*N.y*(N.z/2 + 1));
    
    std::vector<particle> particles = getParticles(p.gets("snapshotFile"), L);
    
    std::cout << "Box dimensions: (" << L.x << ", " << L.y << ", " << L.z << ")" << std::endl;
    std::cout << "Number of particles: " << particles.size() << std::endl;
    
    std::vector<double> frac = {0.25, 0.5, 0.75};
    std::vector<double3> spots;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                double3 pos = {frac[i]*L.x, frac[j]*L.y, frac[k]*L.z};
                spots.push_back(pos);
            }
        }
    }
    
    std::vector<double> kx = fftFreq(N.x, L.x);
    std::vector<double> ky = fftFreq(N.y, L.y);
    std::vector<double> kz = fftFreq(N.z, L.z);
    
    std::random_device seeder;
    std::mt19937_64 gen(seeder());
    std::uniform_real_distribution<double> dist (0.0, 1.0);
    
    particleMeshAcceleration(rho, phi, particles, kx, ky, kz, N, L, p.gets("wisdomFile"), true);
    
    std::ofstream fout(p.gets("outFile"));
    fout.precision(15);
    for (int i = 0; i < 27; ++i) {
        std::cout << "Test #" << i + 1 << ":" << std::endl;
//         double3 pos = {dist(gen)*L.x, dist(gen)*L.y, dist(gen)*L.z};
        
        std::cout << "    Particle-particle calculation..." << std::endl;
        double3 a_pp = getParticleAcceleration(spots[i], particles);
        std::cout << "    Particle-mesh calculation..." << std::endl;
        double3 a_phi = getParticleAcceleration(spots[i], rho, N, L);
        
        double3 diff = a_pp - a_phi;
        
        fout << a_pp.x << " " << a_pp.y << " " << a_pp.z << " " << a_phi.x << " " << a_phi.y << " ";
        fout << a_phi.z << " " << diff.x << " " << diff.y << " " << diff.z << "\n";
    }
    fout.close();
    
    return 0;
}
