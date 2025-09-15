#include "Domain.h"
#include <stdexcept>

MPI_Datatype MPI_DOMAIN;

using namespace std;

ostream& operator<<(ostream& out, const Domain& d)
{
	return out << "Domain[" << d.origin << ", " << d.extent << "]";
}

Domain::Domain()
{
}

Domain::~Domain()
{
}

void Domain::setup(int3 orig, int3 ext)
{
	origin = orig;
	extent = ext;
}

void Domain::BuildMPIDataType()
{
	int block_lengths[50];
	MPI_Aint displacements[50];
	MPI_Aint addresses[50],add_start;
	MPI_Datatype typelist[50];
	Domain temp;
	int count=0;

	typelist[count]=MPI_INT;						// origin.i
	block_lengths[count]=1;
	MPI_Get_address(&temp.origin.i,&addresses[count]);
	count++;
	typelist[count]=MPI_INT;						// origin.j
	block_lengths[count]=1;
	MPI_Get_address(&temp.origin.j,&addresses[count]);
	count++;
	typelist[count]=MPI_INT;						// origin.k
	block_lengths[count]=1;
	MPI_Get_address(&temp.origin.k,&addresses[count]);
	count++;
	typelist[count]=MPI_INT;						// extent.i
	block_lengths[count]=1;
	MPI_Get_address(&temp.extent.i,&addresses[count]);
	count++;
	typelist[count]=MPI_INT;						// extent.j
	block_lengths[count]=1;
	MPI_Get_address(&temp.extent.j,&addresses[count]);
	count++;
	typelist[count]=MPI_INT;						// extent.k
	block_lengths[count]=1;
	MPI_Get_address(&temp.extent.k,&addresses[count]);
	count++;

	MPI_Get_address(&temp,&add_start);
	for (int i=0;i<count;i++) displacements[i]=addresses[i]-add_start;

	MPI_Type_create_struct(count,block_lengths,displacements,typelist,&MPI_DOMAIN);
	MPI_Type_commit(&MPI_DOMAIN);
}

// print to formatted ostream
std::ostream& operator<< (std::ostream& out, const int3& v)
{
	return out << "(" << v.i << "," << v.j << "," << v.k << ")";
}

/*
 * ZFastest specializations
 */
template<>
SubIndex<ZFastest>::SubIndex(const Domain& dom, int array_id)
{
	k = array_id % dom.extent.k + dom.origin.k;
	j = (array_id / dom.extent.k) % dom.extent.j + dom.origin.j;
	i = ((array_id / dom.extent.k) / (dom.extent.j)) % dom.extent.i + dom.origin.i;
}

template<>
int SubIndex<ZFastest>::arrayId(const Domain& dom)
{
	int idx = (k-dom.origin.k) + (j-dom.origin.j)*(dom.extent.k) + (i-dom.origin.i)*(dom.extent.k*dom.extent.j);
	if( !valid(dom) )
	{
		std::stringstream msg;
		msg << "Domain index " << *this << " out of bounds.";
		throw std::runtime_error(msg.str());
	}
	return idx;
}

/*
 * XFastest specializations
 */
template<>
SubIndex<XFastest>::SubIndex(const Domain& dom, int array_id)
{
	i = array_id % dom.extent.i + dom.origin.i;
	j = (array_id / dom.extent.i) % dom.extent.j + dom.origin.j;
	k = ((array_id / dom.extent.i) / (dom.extent.j)) % dom.extent.k + dom.origin.k;
}

template<>
int SubIndex<XFastest>::arrayId(const Domain& dom)
{
	int idx = (i-dom.origin.i) + (j-dom.origin.j)*(dom.extent.i) + (k-dom.origin.k)*(dom.extent.i*dom.extent.j);
	if(idx>dom.extent.size())
	{
		std::stringstream msg;
		msg << "Domain index " << *this << " out of bounds.";
		throw std::runtime_error(msg.str());
	}
	return idx;
}
