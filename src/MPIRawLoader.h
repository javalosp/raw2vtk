#ifndef RAWLOADERMPI_H_
#define RAWLOADERMPI_H_

#include <fstream>
#include <string>
#include <stdexcept>
#include <memory>
#include <algorithm>
#include "MPIDomain.h"

/*
 * Class which loads distinct segments of a RAW voxel
 * image into memory on each process. The class does not
 * actually call any MPI routines as the file is read
 * on each process rather than sending data.
 */
template <typename T, int Padding, IndexScheme S>
class MPIRawLoader : public MPIDomain<T, Padding, S>
{
public:
	MPIRawLoader(std::string fname);
	virtual ~MPIRawLoader();

	void read(size_t header);

private:
	std::string fname;
};

template <typename T, int Padding, IndexScheme S>
MPIRawLoader<T, Padding, S>::MPIRawLoader(std::string fname)
	: fname(fname)
{
}

template <typename T, int Padding, IndexScheme S>
MPIRawLoader<T, Padding, S>::~MPIRawLoader()
{
}

/**
 * @brief Reads a designated segment of a binary .raw file into memory.
 * @details Each MPI process opens the same file but only reads and stores the
 * portion of the data corresponding to its assigned local domain.
 * @param header The size of the file header in bytes to skip before reading voxel data.
 * @throws std::runtime_error if the file cannot be opened.
 */
template <typename T, int Padding, IndexScheme S>
void MPIRawLoader<T, Padding, S>::read(size_t header)
{
	std::ifstream fin(fname.c_str(), std::ios::binary);

	// check file opened okay
	if (!fin.is_open())
		throw std::runtime_error("Cannot open file!");

	// skip header
	fin.seekg(header);

	// read the file, ignoring parts we are not interested in
	// (there might be a quicker way of doing this by reading
	//  whole chunks at a time rather than each number individually)
	T tmp;
	for (int i = 0; i < global.extent.i; i++)
		for (int j = 0; j < global.extent.j; j++)
			for (int k = 0; k < global.extent.k; k++)
			{
				fin.read((char *)&tmp, (std::streamsize)sizeof(T));

				Index idx(i, j, k);

				if (idx.valid(*this))
				{
					MPIDomain<T, Padding, S>::data[idx.arrayId(this->padded)] = tmp;
				}
			}
}

#endif /* RAWLOADERMPI_H_ */
