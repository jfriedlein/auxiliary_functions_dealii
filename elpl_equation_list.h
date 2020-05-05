
#include <deal.II/lac/vector.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include "aux_fncs.h"

using namespace dealii;

template<int dim, typename Number>
class elastoplastic_equations //: public MaterialModel<dim>
 {
 public:
	elastoplastic_equations( const unsigned int plastic_hardening)
	: plastic_hardening(plastic_hardening) {}
	
 private:
	const unsigned int plastic_hardening;
	
	 enum enum_get_elpl
	 {
		 get_alpha = 0,
		 get_R = 1,
		 get_Lambda = 2,
		 get_tangent_plasticAddon = 3,
		 get_dPhi_dgamma = 4
	 };

	 // ToDo: use external enums file
	 enum enum_hard_type
	 {
		 standard_lin_iso = 0,
		 saturated_alpha = 1,
		 saturated_Voce_hard_stress = 2,
		 saturated_VoceLin_hard_stress = 3
	 };
	 
	 // ToDo:
	 // check how to implement/use the parameter like K, beta_inf;
	 // maybe into constructor or use full parameter argument
	 // or use a std::vector that contains all the values and is
	 // accessed only via enumerators to identify which entry corresponds
	 // to which variable
	 const double K = 1;
	 const double beta_inf = 1;
	 const double mu = 1;
	 const double K_exp = 1;
	 
	 // ToDo-optimize: Think about using pointers, etc. to avoid the local copies of 
	 // these variables in the normal code and in this class. Or store the values
	 // only in this class and access them through this class.
	 Number alpha_n_k, gamma, R;
	 double alpha;
	 SymmetricTensor<4,dim> Lambda;
	 SymmetricTensor<4,3> tangent_plasticAddon;
	 SymmetricTensor<2,3,Number> n_n1;
	 Number sigma_dev_t_n1_norm;
	 double dPhi_dgamma;

 public:
   // Accessor functions: \n
   // 1. Update the member variables with the new values from the input arguments \n
   // 2. Call the equation list function and evaluate the desired entry \n
   // 3. Return the updated member variable
	 Number get_alpha_n_k_ep( const double &alpha_n, const Number &gamma_p )
	 {
		alpha = alpha_n;
		gamma = gamma_p;
		elastoplastic_equations_list ( get_alpha );

		return alpha_n_k;
	 };
	 double get_dPhi_dgamma_ep( const double &alpha_n )
	 {
	 	alpha = alpha_n;
	 	elastoplastic_equations_list ( get_dPhi_dgamma );

	 	return dPhi_dgamma;
	 };
	 Number get_hardeningStress_R_ep( const Number &alpha_k )
	 {
	 	alpha_n_k = alpha_k;
	 	elastoplastic_equations_list ( get_R );

	 	return R;
	 };
	 SymmetricTensor<4,dim> get_Lambda_ep( const SymmetricTensor<2,3,Number> &n_n1_in, const Number &gamma_in, const Number &sigma_dev_t_n1_norm_in )
	 {
	 	n_n1 = n_n1_in;
	 	gamma = gamma_in;
	 	sigma_dev_t_n1_norm = sigma_dev_t_n1_norm_in;

	 	elastoplastic_equations_list ( get_Lambda );

	 	return Lambda;
	 };
	 SymmetricTensor<4,dim> get_tangent_plasticAddon_ep( const double &alpha_k, const SymmetricTensor<2,3,Number> &n_n1_in, const Number &gamma_in, const Number &sigma_dev_t_n1_norm_in )
	 {
		alpha = alpha_k;
	 	n_n1 = n_n1_in;
	 	gamma = gamma_in;
	 	sigma_dev_t_n1_norm = sigma_dev_t_n1_norm_in;

	 	elastoplastic_equations_list ( get_tangent_plasticAddon );

	 	return tangent_plasticAddon;
	 };
	 
   // Summary of equations for the different hardening types
	 void elastoplastic_equations_list( enum_get_elpl get_elpl )
	 {
	 	switch ( plastic_hardening )
	 	{
			case standard_lin_iso: { switch ( get_elpl ) {
				//############################################################################################################################################
				case get_alpha: alpha_n_k = alpha + std::sqrt(2./3.) * gamma; break;
				//############################################################################################################################################
				case get_R: R = - K * alpha_n_k; break;
				//############################################################################################################################################
				case get_Lambda:
					 Lambda = 2. * mu / ( 2.*mu + 2./3.*K ) * outer_product(n_n1,n_n1)
							  + ( 2 * mu * gamma ) / sigma_dev_t_n1_norm
							    * ( StandardTensors::I_deviatoric<3>() - outer_product(n_n1,n_n1) ); break;
				//############################################################################################################################################
				case get_tangent_plasticAddon:
					 tangent_plasticAddon = - 4.*mu*mu / -(get_dPhi_dgamma_ep(alpha)) * outer_product(n_n1,n_n1)
											- 4.*mu*mu * gamma / sigma_dev_t_n1_norm
											  * ( StandardTensors::I_deviatoric<3>() - outer_product(n_n1,n_n1) ); break;
				//############################################################################################################################################
				case get_dPhi_dgamma: dPhi_dgamma = - ( 2.*mu + 2./3.*K ); break;
				//############################################################################################################################################
				}
			}
			break;
	 		case saturated_alpha: { switch ( get_elpl ) {
	 			//############################################################################################################################################
				case get_alpha: alpha_n_k = ( alpha + std::sqrt(2./3.) * gamma ) / ( 1. + std::sqrt(2./3.) * K/beta_inf * gamma ); break;
				//############################################################################################################################################
				case get_R: R = - K * alpha_n_k; break;
				//############################################################################################################################################
				case get_Lambda: AssertThrow(false,ExcMessage("elasto-plasticity equation list<< This entry has not yet been implemented.")); break;
				//############################################################################################################################################
				case get_tangent_plasticAddon: ; break;
				//############################################################################################################################################
				case get_dPhi_dgamma: ; break;
				//############################################################################################################################################
				}
	 		}
	 		break;
	 		case saturated_Voce_hard_stress: { switch ( get_elpl ) {
	 			//############################################################################################################################################
				case get_alpha: alpha_n_k = ( alpha + std::sqrt(2./3.) * gamma ) / ( 1. + std::sqrt(2./3.) * K/beta_inf * gamma ); break;
				//############################################################################################################################################
				case get_R: R = - beta_inf * ( 1 - std::exp(-K/beta_inf * alpha_n_k) ); break;
				//############################################################################################################################################
				case get_Lambda: ; break;
				//############################################################################################################################################
				case get_tangent_plasticAddon: ; break;
				//############################################################################################################################################
				case get_dPhi_dgamma: dPhi_dgamma = - (2.*mu + 2./3.*K_exp * std::exp(-K_exp/beta_inf * alpha) ); break;
				//############################################################################################################################################
				}
	 		}
	 		break;
	 		case saturated_VoceLin_hard_stress: { switch ( get_elpl ) {
	 			//############################################################################################################################################
				case get_alpha: alpha_n_k = ( alpha + std::sqrt(2./3.) * gamma ) / ( 1. + std::sqrt(2./3.) * K/beta_inf * gamma ); break;
				//############################################################################################################################################
				case get_R: R = - K * alpha - beta_inf * ( 1 - std::exp(-K_exp * alpha_n_k) ); break;
				//############################################################################################################################################
				case get_Lambda: ; break;
				//############################################################################################################################################
				case get_tangent_plasticAddon: ; break;
				//############################################################################################################################################
				case get_dPhi_dgamma: dPhi_dgamma = - (2.*mu + 2./3.*K + 2./3. * beta_inf * K_exp * std::exp(-K_exp * alpha) ); break;
				//############################################################################################################################################
				}
	 		}
	 		break;
	 		// add case xyz:
	 	}
	 };
 };