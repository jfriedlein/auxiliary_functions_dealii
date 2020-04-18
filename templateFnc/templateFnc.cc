/* ---------------------------------------------------------------------
 *
 * The first coding assignment to get familiar with tensor calculus related
 * deal.II class templates. This includes:
 * 
 * - Tensor<1,dim>
 * - Tensor<2,dim>
 * - Tensor<4,dim>
 * - Vector<double>
 * - FullMatrix<double>
 * 
 * dim is a template variable which allows to vary e.g. between the two- and
 * three dimensional case. As described in the brief repetition of the 
 * essentials in C++, the respective tensor class templates allow different
 * dimensions as input, i.e. dim=1,dim=2,dim=3 for each respective rank
 *
 * ---------------------------------------------------------------------
 */


#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/timer.h>

#include <iostream>
#include <vector>

#include "space.h"

using namespace dealii;



// actual function
template <int dim>
class Solid
{
  public:
	Solid(unsigned int poly_degree);

	unsigned int degree;
	
	void run();
};

template <int dim>
Solid<dim>::Solid( unsigned int poly_degree)
:
degree(poly_degree)
{
}

template <int dim>
void Solid<dim>::run()
{
	std::cout << "RUN" << std::endl;
	
	Tensor<2,dim> F;
	SymmetricTensor<4,dim> C;
	SymmetricTensor<2,dim> hencky_strain,stress;
	
//	SymmetricTensor<4,dim> (*funcPointer) (SymmetricTensor <2,3> &, SymmetricTensor <2,3> &);
//	
//	funcPointer = &ln_space::elastic<dim>;
//	std::cout << "run<< transform elastic" << std::endl;
//	ln_space::transform( F, funcPointer, C);
//	
//	
//	funcPointer = &ln_space::plastic<dim>;
//	std::cout << "run<< transform plastic" << std::endl;
//	ln_space::transform( F, funcPointer, C);
	
	{
		ln_space::ln_space2<dim> ln_space_inst;
		ln_space_inst.pre_ln(F, hencky_strain);
		C = ln_space::elastic<dim> (hencky_strain,stress);
		ln_space_inst.post_ln(stress,C);
	}
	
}



//----------------------------------------------------

int main ()
{
    const int dim=2;
    
    unsigned int degree=2;
    
    Solid<dim> solid_xd(degree);
    solid_xd.run();
}

//----------------------------------------------------