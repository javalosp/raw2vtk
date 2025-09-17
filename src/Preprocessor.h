#ifndef PREPROCESSOR_H_
#define PREPROCESSOR_H_

#include <string>
#include "compiler_opts.h"
#include "MPIDomain.h"
#include "MPIRawLoader.h"

class Preprocessor
{
public:
    Preprocessor();
    virtual ~Preprocessor();

    // Sets up the global domain and decomposes it for each MPI process
    void setupDomain(int3 global_extent);

    // Reads the raw image data from the specified file
    void readRawFile(const std::string &filename, size_t header_size);

    // Writes the material domain to a VTK file set
    void writeVtkFile(const std::string &fname_root);

private:
    template <IndexScheme S>
    void decomposeDomain();

    int mpi_rank;
    int mpi_comm_size;

    Domain local_domain;  // The part of the domain this process owns
    Domain global_domain; // The full simulation domain

    // Data storage for the material types from the RAW file
    MPIDomain<RAWType, 1, IDX_SCHEME> material_data;
};

#endif /* PREPROCESSOR_H_ */