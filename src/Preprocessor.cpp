#include "Preprocessor.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <sys/stat.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkUnsignedShortArray.h>
#include <vtkUnsignedCharArray.h>
#include "MPIDetails.h"

using namespace std;

Preprocessor::Preprocessor()
{
    mpi_rank = MPIDetails::Rank();
    mpi_comm_size = MPIDetails::CommSize();
}

Preprocessor::~Preprocessor() {}

/**
 * @brief Decomposes the global domain among processes.
 * @details This template is specialised for different IndexScheme enums (e.g., ZFastest)
 * to perform a 1D decomposition along the appropriate axis (i.e. ZFastest implies decomposing along X axis,
 * XFastest implies decomposing along Z axis).
 */
template <>
void Preprocessor::decomposeDomain<ZFastest>()
{
    local_domain = global_domain;
    int block_size = global_domain.extent.i / mpi_comm_size;
    local_domain.origin.i = mpi_rank * block_size;

    // This ensures the last process gets all remaining voxels
    local_domain.extent.i = (mpi_rank == mpi_comm_size - 1) ? (global_domain.extent.i - local_domain.origin.i) : block_size;
}

template <>
void Preprocessor::decomposeDomain<XFastest>()
{
    local_domain = global_domain;
    int block_size = global_domain.extent.k / mpi_comm_size;
    local_domain.origin.k = mpi_rank * block_size;
    local_domain.extent.k = (mpi_rank == mpi_comm_size - 1) ? (global_domain.extent.k - local_domain.origin.k) : block_size;
}

/**
 * @brief Sets up the global simulation domain and decomposes it across MPI processes.
 * @param gextent The dimensions (i, j, k) of the entire dataset.
 */
void Preprocessor::setupDomain(int3 gextent)
{
    global_domain.origin = int3();
    global_domain.extent = gextent;

    decomposeDomain<IDX_SCHEME>();

    MPIDomain<double, 0, IDX_SCHEME>::SetGlobal(int3(), gextent);
    MPISubIndex<IDX_SCHEME>::Init(local_domain, mpi_rank, mpi_comm_size);

    material_data.setup(local_domain.origin, local_domain.extent);

    if (mpi_rank == 0)
    {
        std::cout << "Global domain setup complete: " << global_domain.extent << std::endl;
    }
}

/**
 * @brief Manages the reading of the raw image file into the distributed domain.
 * @details Verifies that the file size matches the expected domain size and then
 * uses MPIRawLoader to perform the parallel read.
 * @param filename The path to the .raw input file.
 * @param header_size The size of the file header in bytes.
 * @throws std::runtime_error if the file size does not match the domain dimensions.
 */
void Preprocessor::readRawFile(const std::string &filename, size_t header_size)
{
    if (mpi_rank == 0)
    {
        std::cout << "Reading RAW file: " << filename << std::endl;
    }

    struct stat filestatus;
    if (stat(filename.c_str(), &filestatus) != 0)
    {
        throw std::runtime_error("Cannot get file status for " + filename);
    }

    if ((filestatus.st_size - header_size) != global_domain.extent.size() * sizeof(RAWType))
    {
        std::stringstream msg;
        msg << "File size does not match specified domain dimensions." << std::endl;
        msg << "\tFile size on disk: " << filestatus.st_size << " bytes." << std::endl;
        msg << "\tExpected data size: " << global_domain.extent.size() * sizeof(RAWType) << " bytes.";
        throw std::runtime_error(msg.str());
    }

    MPIRawLoader<RAWType, 1, IDX_SCHEME> reader(filename);
    reader.setup(local_domain.origin, local_domain.extent);
    reader.read(header_size);
    material_data.take(reader.getData());

    if (mpi_rank == 0)
    {
        std::cout << "RAW file reading complete." << std::endl;
    }
}

/**
 * @brief Writes the data to a set of VTK files.
 * @details The root process writes a master .pvti file that references individual .vti
 * part files written by each process.
 * @param fname_root The base filename for the output files (e.g., "./output/material").
 */
void Preprocessor::writeVtkFile(const std::string &fname_root)
{
    // Write the master .pvti file on the root process
    if (mpi_rank == 0)
    {
        std::stringstream pvti_fname;
        pvti_fname << fname_root << ".pvti";
        std::ofstream fout(pvti_fname.str());

        fout << "<?xml version=\"1.0\"?>" << std::endl;
        fout << "<VTKFile type=\"PImageData\" version=\"0.1\">" << std::endl;
        fout << "\t<PImageData WholeExtent=\"0 " << global_domain.extent.i - 1
             << " 0 " << global_domain.extent.j - 1 << " 0 "
             << global_domain.extent.k - 1 << "\" ";
        fout << "GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
        fout << "\t\t<PPointData Scalars=\"MaterialType\">" << std::endl;

#if DATA_TYPE == 16
        fout << "\t\t\t<PDataArray type=\"UInt16\" Name=\"MaterialType\"/>" << std::endl;
#elif DATA_TYPE == 8
        fout << "\t\t\t<PDataArray type=\"UInt8\" Name=\"MaterialType\"/>" << std::endl;
#endif

        fout << "\t\t</PPointData>" << std::endl;

        // Write the references to the part files with their corresponding extent
        for (int proc = 0; proc < mpi_comm_size; ++proc)
        {
            std::stringstream piece_fname;
            std::string root_basename = fname_root.substr(fname_root.find_last_of("/\\") + 1);
            piece_fname << root_basename << "_" << proc << ".vti";

            const Domain &piece_dom = MPISubIndex<IDX_SCHEME>::all_local_domains[proc];
            int max_i = piece_dom.origin.i + piece_dom.extent.i - 1;

            // Add the one-voxel overlap for all pieces except the very last one
            if (proc != mpi_comm_size - 1)
            {
                max_i++;
            }

            fout << "\t\t<Piece Extent=\"";
            fout << piece_dom.origin.i << " " << max_i << " ";
            fout << piece_dom.origin.j << " " << piece_dom.origin.j + piece_dom.extent.j - 1 << " ";
            fout << piece_dom.origin.k << " " << piece_dom.origin.k + piece_dom.extent.k - 1 << "\" ";
            fout << "Source=\"" << piece_fname.str() << "\"/>" << std::endl;
        }
        fout << "\t</PImageData>" << std::endl;
        fout << "</VTKFile>" << std::endl;
        fout.close();
    }

    // Ensure all processes wait for rank 0 to finish writing the master file
    MPI_Barrier(MPI_COMM_WORLD);

    // Each process writes its own .vti part file using the "safe" copy method
    std::stringstream vti_fname;
    vti_fname << fname_root << "_" << mpi_rank << ".vti";

    int off_pos = (mpi_rank == mpi_comm_size - 1) ? 0 : 1;

    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetExtent(local_domain.origin.i, local_domain.origin.i + local_domain.extent.i - 1 + off_pos,
                         local_domain.origin.j, local_domain.origin.j + local_domain.extent.j - 1,
                         local_domain.origin.k, local_domain.origin.k + local_domain.extent.k - 1);

    size_t num_voxels_to_write = (size_t)(local_domain.extent.i + off_pos) * local_domain.extent.j * local_domain.extent.k;

    // Create the appropriate VTK array and allocate memory inside it
#if DATA_TYPE == 16
    vtkSmartPointer<vtkUnsignedShortArray> type_arr = vtkSmartPointer<vtkUnsignedShortArray>::New();
#elif DATA_TYPE == 8
    vtkSmartPointer<vtkUnsignedCharArray> type_arr = vtkSmartPointer<vtkUnsignedCharArray>::New();
#endif
    type_arr->SetName("MaterialType");
    type_arr->SetNumberOfComponents(1);
    type_arr->SetNumberOfValues(num_voxels_to_write);

    // Loop through all voxels to be written (including the overlap)
    // and copy the data from your simulation buffer to the VTK buffer.
    size_t count = 0;
    for (int k = local_domain.origin.k; k < (local_domain.origin.k + local_domain.extent.k); ++k)
    {
        for (int j = local_domain.origin.j; j < (local_domain.origin.j + local_domain.extent.j); ++j)
        {
            for (int i = local_domain.origin.i; i < (local_domain.origin.i + local_domain.extent.i + off_pos); ++i)
            {
                Index idx(i, j, k);
                type_arr->SetValue(count, material_data[idx]);
                count++;
            }
        }
    }

    imageData->GetPointData()->AddArray(type_arr);

    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(vti_fname.str().c_str());
    writer->SetInputConnection(imageData->GetProducerPort());
    writer->Write();
}
