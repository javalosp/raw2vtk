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

Preprocessor::Preprocessor()
{
    mpi_rank = MPIDetails::Rank();
    mpi_comm_size = MPIDetails::CommSize();
}

Preprocessor::~Preprocessor() {}

template <>
void Preprocessor::decomposeDomain<ZFastest>()
{
    local_domain = global_domain;
    local_domain.origin.i = mpi_rank * (global_domain.extent.i / mpi_comm_size);
    local_domain.extent.i = global_domain.extent.i / mpi_comm_size;
    if (mpi_rank == mpi_comm_size - 1)
    {
        local_domain.extent.i = global_domain.extent.i - (local_domain.origin.i - global_domain.origin.i);
    }
}

template <>
void Preprocessor::decomposeDomain<XFastest>()
{
    local_domain = global_domain;
    local_domain.origin.k = mpi_rank * (global_domain.extent.k / mpi_comm_size);
    local_domain.extent.k = global_domain.extent.k / mpi_comm_size;
    if (mpi_rank == mpi_comm_size - 1)
    {
        local_domain.extent.k = global_domain.extent.k - (local_domain.origin.k - global_domain.origin.k);
    }
}

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

void Preprocessor::writeVtkFile(const std::string &fname_root)
{
    // Write the master .pvti file from the root process
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

        // --- START FIX: Use the new DATA_TYPE macro ---
#if DATA_TYPE == 16
        fout << "\t\t\t<PDataArray type=\"UInt16\" Name=\"MaterialType\"/>" << std::endl;
#elif DATA_TYPE == 8
        fout << "\t\t\t<PDataArray type=\"UInt8\" Name=\"MaterialType\"/>" << std::endl;
#endif
        // --- END FIX ---

        fout << "\t\t</PPointData>" << std::endl;

        for (int proc = 0; proc < mpi_comm_size; ++proc)
        {
            std::stringstream vti_fname;
            vti_fname << fname_root << "_" << proc << ".vti";
            fout << "\t\t<Piece Extent=\"";
            fout << MPISubIndex<IDX_SCHEME>::all_local_domains[proc].origin.i << " "
                 << MPISubIndex<IDX_SCHEME>::all_local_domains[proc].origin.i + MPISubIndex<IDX_SCHEME>::all_local_domains[proc].extent.i - 1 << " ";
            fout << MPISubIndex<IDX_SCHEME>::all_local_domains[proc].origin.j << " "
                 << MPISubIndex<IDX_SCHEME>::all_local_domains[proc].origin.j + MPISubIndex<IDX_SCHEME>::all_local_domains[proc].extent.j - 1 << " ";
            fout << MPISubIndex<IDX_SCHEME>::all_local_domains[proc].origin.k << " "
                 << MPISubIndex<IDX_SCHEME>::all_local_domains[proc].origin.k + MPISubIndex<IDX_SCHEME>::all_local_domains[proc].extent.k - 1 << "\" ";
            fout << "Source=\"" << vti_fname.str().substr(fname_root.find_last_of("/\\") + 1) << "\"/>" << std::endl;
        }
        fout << "\t</PImageData>" << std::endl;
        fout << "</VTKFile>" << std::endl;
        fout.close();
        std::cout << "Wrote master VTK file: " << pvti_fname.str() << std::endl;
    }

    // Each process writes its own .vti piece file
    std::stringstream vti_fname;
    vti_fname << fname_root << "_" << mpi_rank << ".vti";

    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetExtent(
        local_domain.origin.i, local_domain.origin.i + local_domain.extent.i - 1,
        local_domain.origin.j, local_domain.origin.j + local_domain.extent.j - 1,
        local_domain.origin.k, local_domain.origin.k + local_domain.extent.k - 1);

    // --- START FIX: Use the new DATA_TYPE macro to select the correct VTK Array type ---
#if DATA_TYPE == 16
    vtkSmartPointer<vtkUnsignedShortArray> type_arr = vtkSmartPointer<vtkUnsignedShortArray>::New();
#elif DATA_TYPE == 8
    vtkSmartPointer<vtkUnsignedCharArray> type_arr = vtkSmartPointer<vtkUnsignedCharArray>::New();
#endif
    // --- END FIX ---

    type_arr->SetName("MaterialType");
    type_arr->SetNumberOfComponents(1);

    type_arr->SetArray(material_data.getData().get() + material_data.pad_size, local_domain.extent.size(), 1);

    imageData->GetPointData()->AddArray(type_arr);

    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(vti_fname.str().c_str());
    writer->SetInputConnection(imageData->GetProducerPort());
    writer->Write();
}