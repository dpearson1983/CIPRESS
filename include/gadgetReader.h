#ifndef _GADGETREADER_H_
#define _GADGETREADER_H_

#include <vector>
#include <string>

// Convenience structure for storing vector data
struct float3{
    float x, y, z;
};

// Data structure for GADGET2's custom header in the snapshot files. "unused" is empty space left
// by GADGET2's authors for possible adding additional items later.
struct gadget_header{
    unsigned int N_p[6];
    double M_p[6];
    double time;
    double redshift;
    int flagSfr;
    int flagFeedback;
    int N_tot[6];
    int flagCooling;
    int numFiles;
    double boxSize;
    double Omega_0;
    double Omega_L;
    double h_0;
    int flagAge;
    int flagMetals;
    unsigned int N_totHW[6];
    int flagEntrICs;
    char unused[256 - 6*4 - 6*8 - 2*8 - 2*4 - 6*4 - 2*4 - 4*8 - 2*4 - 6*4 - 4];
};
    

void read_gadget_snapshot(std::string snapshot, gadget_header &header, 
                          std::vector<std::vector<float3>> &pos,
                          std::vector<std::vector<float3>> &vel);

void read_gadget_snapshot(std::string snapshot, gadget_header &header, 
                          std::vector<std::vector<float3>> &pos,
                          std::vector<std::vector<float3>> &vel, std::vector<std::vector<int>> &ids);

void print_gadget_header(gadget_header &header);

#endif
