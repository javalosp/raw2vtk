#include <iostream>
#include <stdexcept>
#include <mpi.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "compiler_opts.h"
#include "Domain.h"
#include "Preprocessor.h"

namespace opts = boost::program_options;

/**
 * @brief Main entry point for the RAW to VTK preprocessing application.
 * @details This function executes the preprocessing workflow. It initialises MPI,
 * parses command-line arguments for domain dimensions and file paths, sets up the
 * distributed computational domain, reads the source .raw file, and writes the
 * output to a parallel VTK file format.
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line arguments.
 * @return Returns 0 on successful execution, 1 on error.
 */
int main(int argc, char *argv[])
{
    try
    {
        MPI_Init(&argc, &argv);

        int mpi_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

        // Command line arguments
        opts::options_description cmd_opts("Usage");
        cmd_opts.add_options()("help,h", "Print this help message")("raw-file", opts::value<std::string>()->required(), "Input RAW file specifying the domain.")("x-ext", opts::value<int>()->required(), "The x extent (width) of the domain.")("y-ext", opts::value<int>()->required(), "The y extent (height) of the domain.")("z-ext", opts::value<int>()->required(), "The z extent (depth) of the domain.")("header-size", opts::value<size_t>()->default_value(0), "RAW file header size in bytes.")("output-dir", opts::value<std::string>()->default_value("./output"), "The output directory for VTK files.");

        opts::variables_map vm;
        try
        {
            opts::store(opts::parse_command_line(argc, argv, cmd_opts), vm);

            if (vm.count("help"))
            {
                if (mpi_rank == 0)
                {
                    std::cout << "RAW to VTK Preprocessor" << std::endl;
                    std::cout << cmd_opts << std::endl;
                }
                MPI_Finalize();
                return 0;
            }

            opts::notify(vm);
        }
        catch (const opts::error &e)
        {
            if (mpi_rank == 0)
            {
                std::cerr << "Error parsing command line options: " << e.what() << std::endl;
                std::cerr << cmd_opts << std::endl;
            }
            MPI_Finalize();
            return 1;
        }

        // Setup MPI and Domain data types
        Domain::BuildMPIDataType();

        // Create and run the preprocessor
        Preprocessor preprocessor;

        // int3 global_extent(vm["x-ext"].as<int>(), vm["y-ext"].as<int>(), vm["z-ext"].as<int>());
        // Map arguments to the code's (i, j, k) = (Z, Y, X) internal indexing
        int3 global_extent(vm["z-ext"].as<int>(), vm["y-ext"].as<int>(), vm["x-ext"].as<int>());

        preprocessor.setupDomain(global_extent);
        preprocessor.readRawFile(vm["raw-file"].as<std::string>(), vm["header-size"].as<size_t>());

        // Ensure the output directory exists
        std::string out_dir = vm["output-dir"].as<std::string>();

        // Create the directory if it doesn't exist (only on rank 0 to avoid race conditions)
        if (mpi_rank == 0)
        {
            boost::filesystem::create_directories(out_dir);
        }
        // Ensure all processes wait until the directory is created before proceeding
        MPI_Barrier(MPI_COMM_WORLD);

        // Write output files
        preprocessor.writeVtkFile(out_dir + "/material_domain");

        MPI_Finalize();
    }
    catch (const std::exception &e)
    {
        std::cerr << "An unhandled exception occurred: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    catch (...)
    {
        std::cerr << "An unknown error occurred." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    return 0;
}