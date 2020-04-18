#ifndef ln_space_H
#define ln_space_H

// @section includes Include Files
// The data type SymmetricTensor and some related operations, such as trace, symmetrize, deviator, ... for tensor calculus
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/exceptions.h>

// @todo Check whether the following three headers are needed at all
#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

namespace ln_space
{
//	template<int dim>
//	void transform (  /*input->*/ Tensor<2, dim> &F,
//									SymmetricTensor<4,3> (*func_pointer)(SymmetricTensor <2,3> &, SymmetricTensor <2,3> & ),
//			 	 	  /*output->*/ SymmetricTensor<4,dim> &C )
//	{
//		std::cout << "ln_space::transform<< running" << std::endl;
//
//		 SymmetricTensor<2, dim> hencky_strain, T_stress;
//		 
//		 C = func_pointer( /*input->*/ hencky_strain, T_stress );
//	}

	// material model interface: input-> strain, phi, history_n
	//							 output-> stress, vM_stress, history_tmp, C + C_ij + Lambda + GG_mode_requested + hardening type + material_id + eval_QP
	

// actual function
template <int dim>
class ln_space2
{
  public:
	ln_space2();

	Vector<double> abc;
	
	void pre_ln ( 	/*input->*/ Tensor<2, dim> &F, 
					/*output->*/ SymmetricTensor<2,dim> &hencky_strain  );

	void post_ln ( /*output->*/ SymmetricTensor<2,dim> &second_piola_stress_S, SymmetricTensor<4,dim> &elasto_plastic_tangent	);

};

template <int dim>
ln_space2<dim>::ln_space2( )
:
abc (3)
{
}


template<>
void ln_space2<3>::pre_ln ( 	/*input->*/ Tensor<2, 3> &F, 
				/*output->*/ SymmetricTensor<2,3> &hencky_strain  ) 
{
	std::cout << "ln_space::pre_ln<< running" << std::endl;

	abc(0) = 33;
	abc(1) = 34;
}


template<>
void ln_space2<2>::pre_ln ( 	/*input->*/ Tensor<2, 2> &F, 
				/*output->*/ SymmetricTensor<2,2> &hencky_strain  ) 
{
	std::cout << "ln_space::pre_ln<< running" << std::endl;

	abc(0) = 32;
	abc(1) = 322;
}

template<int dim>
void ln_space2<dim>::post_ln ( /*output->*/ SymmetricTensor<2,dim> &second_piola_stress_S, SymmetricTensor<4,dim> &elasto_plastic_tangent	)
{
	std::cout << "ln_space::post_ln<< running" << std::endl;

	std::cout << "abc=" << abc << std::endl;
}

	template<int dim>
	SymmetricTensor<4,2> elastic ( SymmetricTensor <2,2> &strain, SymmetricTensor <2,2> &stress  )
	{
		std::cout << "ln_space::elastic<< running" << std::endl;
		return outer_product(strain,stress);
	}
	
	template<int dim>
	SymmetricTensor<4,3> plastic ( SymmetricTensor <2,3> &strain, SymmetricTensor <2,3> &stress  )
	{
		std::cout << "ln_space::plastic<< running" << std::endl;
		return outer_product(strain,stress);
	}
}


#endif // ln_space_H