#include "../include/gadgetReader.h"
#include <iostream>
#include <fstream>
#include <sstream>

std::vector<std::string> part_names = {"Gas", "Halo", "Disk", "Bulge", "Star", "Boundary"};

void read_gadget_snapshot(std::string snapshot, gadget_header &header, 
                          std::vector<std::vector<float3>> &pos,
                          std::vector<std::vector<float3>> &vel) {
    std::ifstream fin(snapshot, std::ios::in|std::ios::binary);
    
    int blockSize;
    
    // Read all the crap in the header of the file
    fin.read((char *) &blockSize, sizeof(int));
    fin.read((char *) &header, sizeof(header));
    fin.read((char *) &blockSize, sizeof(int));
    
    // Get the size of the position block
    int posBlockSize;
    fin.read((char *) &posBlockSize, sizeof(int));
    
    int theoPosBlockSize = 0;
    for (int i = 0; i < 6; ++i)
        theoPosBlockSize += header.N_p[i]*sizeof(float3);
    
    if (posBlockSize != theoPosBlockSize) {
        print_gadget_header(header);
        std::stringstream errMessage;
        errMessage << "The reported size of the position block does not appear to be correct.";
        throw std::runtime_error(errMessage.str());
    }
    
    // Read in particle positions
    for (int i = 0; i < 6; ++i) {
        if (header.N_p[i] > 0) {
            for (unsigned int j = 0; j < header.N_p[i]; ++j) {
                float3 part;
                fin.read((char *) &part, sizeof(float3));
                pos[i].push_back(part);
            }
        }
    }
    
    // Verify end of position block
    fin.read((char *) &blockSize, sizeof(int));
    
    // Get the size of the velocity block
    int velBlockSize;
    fin.read((char *) &velBlockSize, sizeof(int));
    
    int theoVelBlockSize = 0;
    for (int i = 0; i < 6; ++i)
        theoVelBlockSize += header.N_p[i]*sizeof(float3);
    
    if (velBlockSize != theoVelBlockSize) {
        print_gadget_header(header);
        std::stringstream errMessage;
        errMessage << "The reported size of the velocity block does not appear to be correct.";
        throw std::runtime_error(errMessage.str());
    }
    
    // Read in the particle velocities
    for (int i = 0; i < 6; ++i) {
        if (header.N_p[i] > 0) {
            for (unsigned int j = 0; j < header.N_p[i]; ++j) {
                float3 part;
                fin.read((char *) &part, sizeof(float3));
                vel[i].push_back(part);
            }
        }
    }
    
    fin.close();
    print_gadget_header(header);
}

void read_gadget_snapshot(std::string snapshot, gadget_header &header, 
                          std::vector<std::vector<float3>> &pos,
                          std::vector<std::vector<float3>> &vel, std::vector<std::vector<int>> &ids) {
    std::ifstream fin(snapshot, std::ios::in|std::ios::binary);
    
    int blockSize;
    
    // Read all the crap in the header of the file
    fin.read((char *) &blockSize, sizeof(int));
    fin.read((char *) &header, sizeof(header));
    fin.read((char *) &blockSize, sizeof(int));
    
    // Get the size of the position block
    int posBlockSize;
    fin.read((char *) &posBlockSize, sizeof(int));
    
    int theoPosBlockSize = 0;
    for (int i = 0; i < 6; ++i)
        theoPosBlockSize += header.N_p[i]*sizeof(float3);
    
    if (posBlockSize != theoPosBlockSize) {
        print_gadget_header(header);
        std::stringstream errMessage;
        errMessage << "The reported size of the position block does not appear to be correct.";
        throw std::runtime_error(errMessage.str());
    }
    
    // Read in particle positions
    for (int i = 0; i < 6; ++i) {
        if (header.N_p[i] > 0) {
            for (unsigned int j = 0; j < header.N_p[i]; ++j) {
                float3 part;
                fin.read((char *) &part, sizeof(float3));
                pos[i].push_back(part);
            }
        }
    }
    
    // Verify end of position block
    fin.read((char *) &blockSize, sizeof(int));
    
    // Get the size of the velocity block
    int velBlockSize;
    fin.read((char *) &velBlockSize, sizeof(int));
    
    int theoVelBlockSize = 0;
    for (int i = 0; i < 6; ++i)
        theoVelBlockSize += header.N_p[i]*sizeof(float3);
    
    if (velBlockSize != theoVelBlockSize) {
        print_gadget_header(header);
        std::stringstream errMessage;
        errMessage << "The reported size of the velocity block does not appear to be correct.";
        throw std::runtime_error(errMessage.str());
    }
    
    // Read in the particle velocities
    for (int i = 0; i < 6; ++i) {
        if (header.N_p[i] > 0) {
            for (unsigned int j = 0; j < header.N_p[i]; ++j) {
                float3 part;
                fin.read((char *) &part, sizeof(float3));
                vel[i].push_back(part);
            }
        }
    }
    
    // Verify the end of the velocity block
    fin.read((char *) &blockSize, sizeof(int));
    
    // Get the size of the particle ID block
    int idBlockSize;
    fin.read((char *) &idBlockSize, sizeof(int));
    
    int theoIDBlockSize = 0;
    for (int i = 0; i < 6; ++i)
        theoIDBlockSize += header.N_p[i]*sizeof(int);
    
    if (idBlockSize != theoIDBlockSize) {
        print_gadget_header(header);
        std::stringstream errMessage;
        errMessage << "The reported size of the particle ID block does not appear to be correct.";
        throw std::runtime_error(errMessage.str());
    }
    
    // Read in the particle IDs
    for (int i = 0; i < 6; ++i) {
        if (header.N_p[i] > 0) {
            for (unsigned int j = 0; j < header.N_p[i]; ++j) {
                int ID;
                fin.read((char *) &ID, sizeof(int));
                ids[i].push_back(ID);
            }
        }
    }
    
    fin.close();
    print_gadget_header(header);
}

void print_gadget_header(gadget_header &header) {
    std::cout << "Snapshot Details:\n";
    std::cout << "      a = " << header.time << " (z = " << header.redshift << ")\n";
    std::cout << "      L = " << header.boxSize << " h^-1 kpc\n";
    std::cout << "Omega_0 = " << header.Omega_0 << "\n";
    std::cout << "Omega_L = " << header.Omega_L << "\n";
    std::cout << "    h_0 = " << header.h_0 << "\n\n";
    std::cout << "Particle Numbers:\n";
    for (int i = 0; i < 6; ++i)
        std::cout << "    Number of " << part_names[i] << " particles: " << header.N_p[i] << " of " 
                  << header.N_tot[i] << "\n";
    std::cout << "\nParticle Masses:\n";
    for (int i = 0; i < 6; ++i)
        std::cout << "    M_" << part_names[i] << " = " << header.M_p[i] << "\n";
}
