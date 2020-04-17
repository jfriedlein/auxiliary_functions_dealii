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

using namespace dealii;


template<int dim>
SymmetricTensor<4,dim> outer_product_sym( const SymmetricTensor<2,dim> &A, const SymmetricTensor<2,dim> &B )
{
	SymmetricTensor<4,dim> D;
	// Special nested for-loop to access only non-symmetric entries of 4th order sym. tensor
	// ToDo: still not optimal element 1112 and 1211 are both accessed
	for ( unsigned int i=0; i<dim; ++i )
		for ( unsigned int j=i; j<dim; ++j )
			for ( unsigned int k=i; k<dim; ++k )
				for ( unsigned int l=k; l<dim; ++l )
				{
					double tmp = A[i][j] * B[k][l] + B[i][j] * A[k][l];
					D[i][j][k][l] = tmp;
					D[k][l][i][j] = tmp;
				}
	return D;
}
template<int dim>
SymmetricTensor<4,dim> outer_product_sym( const SymmetricTensor<2,dim> &A )
{
	SymmetricTensor<4,dim> D;
	// Special nested for-loop to access only non-symmetric entries of 4th order sym. tensor
	// ToDo: still not optimal element 1112 and 1211 are both accessed
	for ( unsigned int i=0; i<dim; ++i )
    	for ( unsigned int j=i; j<dim; ++j )
        	for ( unsigned int k=i; k<dim; ++k )
            	for ( unsigned int l=k; l<dim; ++l )
            	{
            		double tmp = A[i][j] * A[k][l];
            		D[i][j][k][l] = tmp;
            		D[k][l][i][j] = tmp;
            	}
	return D;
}


//----------------------------------------------------

int main ()
{
    /* For now the three dimensional case is considered.
     * This information is essential for deal.II since
     * is suited for arbitrary dimensions due to its
     * template character. Use this variable for all
     * deal.II templates that will be created in the
     * sequent
     * 
     * It is always helpful to read the manual and see
     * how functions are implemented, i.e. what is the
     * return value, how is the function called or what
     * is the input.
     * 
     * Further some functions may be declared 
     * "DEPRECATED" which means they are still usable
     * but will be removed in future releases of the 
     * library -> not recommended to use those
     */
    
    
    const int dim=3;
    //------------------------------------------------
    /* START IMPLEMENTATION HERE                    */
    //------------------------------------------------
    
    
    //------------------------------------------------
    //          EX - 1
    /* Create two tensors of rank one and name them
     * u and v respectively and print them to the screen.
     * Therefore consider the available documentation
     * and manual on the deal.II website
     */
    TimerOutput timer (std::cout, TimerOutput::summary,
                       TimerOutput::cpu_times);
    
    Tensor<1,dim> u;
    Tensor<1,dim> v;
    Tensor<1,dim> w;
    
    //BEGIN - INSERT YOUR CODE HERE
    u[0]=1; u[1]=2; u[2]=3;        
    v[0]=4; v[1]=5; v[2]=6;
    w[0]=7; w[1]=8; w[2]=9;
    //END - INSERT YOUR CODE HERE
    
    std::vector< SymmetricTensor<2,dim> > Ma (3);
    Ma[0] = symmetrize(outer_product(u,u));
    Ma[1] = symmetrize(outer_product(v,v));
    Ma[2] = symmetrize(outer_product(w,w));

    std::cout << Ma[0] << std::endl;
    std::cout << Ma[1] << std::endl;
    std::cout << Ma[2] << std::endl;


    // in RELEASE mode
    // ToDo: for reliable numbers we need to average
    
    // @section Sym2 Outer product of two identical 2nd order symmetric tensors
    
    // This takes 0.0059 s
    SymmetricTensor<4,dim> Ma_x_Ma;
    timer.enter_subsection ("Ma_x_Ma");
    for ( unsigned int i=0; i<10000; i++)
    	Ma_x_Ma = outer_product(Ma[0],Ma[0]);
    timer.leave_subsection ("Ma_x_Ma");

    // This takes 0.0225 s
    SymmetricTensor<4,dim> Ma_x_Ma_for;
    timer.enter_subsection ("Ma_x_Ma-for");
    for ( unsigned int i=0; i<10000; i++)
    {
    	for ( unsigned int i=0; i<dim; ++i )
        	for ( unsigned int j=0; j<dim; ++j )
            	for ( unsigned int k=0; k<dim; ++k )
                	for ( unsigned int l=0; l<dim; ++l )
                		Ma_x_Ma_for[i][j][k][l] = Ma[0][i][j] * Ma[0][k][l];
    }
    timer.leave_subsection ("Ma_x_Ma-for");
    
    // Trying to exploit symmetry properties:
    // "This entails certain symmetry properties on the elements in their 4-dimensional index space,
    // in particular that Cijkl=Cjikl=Cijlk. However, it does not imply the relation Cijkl=Cklij"
    // This takes 0.00485 s
    SymmetricTensor<4,dim> Ma_x_Ma_for2;
    timer.enter_subsection ("Ma_x_Ma-for2");
    for ( unsigned int i=0; i<10000; i++)
    {
    	// Special nested for-loop to access only non-symmetric entries of 4th order sym. tensor
    	for ( unsigned int i=0; i<dim; ++i )
        	for ( unsigned int j=i; j<dim; ++j )
            	for ( unsigned int k=0; k<dim; ++k )
                	for ( unsigned int l=k; l<dim; ++l )
                		Ma_x_Ma_for2[i][j][k][l] = Ma[0][i][j] * Ma[0][k][l];
    }
    timer.leave_subsection ("Ma_x_Ma-for2");
    
    // BUT: Due to our (Ma x Ma) setup we even have major symmetry
    // This takes 0.00391 s
    SymmetricTensor<4,dim> Ma_x_Ma_for3;
    timer.enter_subsection ("Ma_x_Ma-for3");
    for ( unsigned int i=0; i<10000; i++)
    {
    	// Special nested for-loop to access only non-symmetric entries of 4th order sym. tensor
    	// ToDo: still not optimal element 1112 and 1211 are both accessed
    	for ( unsigned int i=0; i<dim; ++i )
        	for ( unsigned int j=i; j<dim; ++j )
            	for ( unsigned int k=i; k<dim; ++k )
                	for ( unsigned int l=k; l<dim; ++l )
                	{
                		double tmp = Ma[0][i][j] * Ma[0][k][l];
                		Ma_x_Ma_for3[i][j][k][l] = tmp;
                		Ma_x_Ma_for3[k][l][i][j] = tmp;
                	}
    }
    timer.leave_subsection ("Ma_x_Ma-for3");
    
    SymmetricTensor<4,dim> Ma_x_Ma_for4;
    timer.enter_subsection ("Ma_x_Ma-for4");
    for ( unsigned int i=0; i<10000; i++)
    	Ma_x_Ma_for4 = outer_product_sym(Ma[0]);
    timer.leave_subsection ("Ma_x_Ma-for4");
    
    double error_for = 0.;
    double error_for2 = 0.;
    double error_for3 = 0.;
    double error_for4 = 0.;
	for ( unsigned int i=0; i<dim; ++i )
    	for ( unsigned int j=0; j<dim; ++j )
        	for ( unsigned int k=0; k<dim; ++k )
            	for ( unsigned int l=0; l<dim; ++l )
            	{
            		error_for += std::abs(Ma_x_Ma[i][j][k][l]-Ma_x_Ma_for[i][j][k][l]);
            		error_for2 += std::abs(Ma_x_Ma[i][j][k][l]-Ma_x_Ma_for2[i][j][k][l]);
            		error_for3 += std::abs(Ma_x_Ma[i][j][k][l]-Ma_x_Ma_for3[i][j][k][l]);
            		error_for4 += std::abs(Ma_x_Ma[i][j][k][l]-Ma_x_Ma_for4[i][j][k][l]);
            	}
	
	std::cout << std::endl;
    std::cout << "Ma_x_Ma=     " << Ma_x_Ma << std::endl;
    std::cout << "Ma_x_Ma_for= " << Ma_x_Ma_for << " with error=" << error_for << std::endl;
    std::cout << "Ma_x_Ma_for2=" << Ma_x_Ma_for2 << " with error=" << error_for2 << std::endl;
    std::cout << "Ma_x_Ma_for3=" << Ma_x_Ma_for3 << " with error=" << error_for3 << std::endl;
    std::cout << "Ma_x_Ma_for4=" << Ma_x_Ma_for4 << " with error=" << error_for4 << std::endl;

    
    // @section Sym2 Outer product of two identical 2nd order symmetric tensors
    SymmetricTensor<4,dim> D;
    timer.enter_subsection ("D");
    for ( unsigned int i=0; i<10000; i++)
    	D = outer_product(Ma[0],Ma[1]) + outer_product(Ma[1],Ma[0]);
    timer.leave_subsection ("D");
    
    SymmetricTensor<4,dim> D_for1;
    timer.enter_subsection ("D_for1");
    for ( unsigned int i=0; i<10000; i++)
    {
    	// Special nested for-loop to access only non-symmetric entries of 4th order sym. tensor
    	// ToDo: still not optimal element 1112 and 1211 are both accessed
    	for ( unsigned int i=0; i<dim; ++i )
        	for ( unsigned int j=i; j<dim; ++j )
            	for ( unsigned int k=i; k<dim; ++k )
                	for ( unsigned int l=k; l<dim; ++l )
                	{
                		double tmp = Ma[0][i][j] * Ma[1][k][l] + Ma[1][i][j] * Ma[0][k][l];
                		D_for1[i][j][k][l] = tmp;
                		D_for1[k][l][i][j] = tmp;
                	}
    }
    timer.leave_subsection ("D_for1");

    SymmetricTensor<4,dim> D_for2;
    timer.enter_subsection ("D_for2");
    for ( unsigned int i=0; i<10000; i++)
    	D_for2 = outer_product_sym(Ma[0], Ma[1]);	// the same as (Ma[1],Ma[0])
    timer.leave_subsection ("D_for2");
    
    double error_dfor = 0.;
    double error_dfor2 = 0.;
//    double error_dfor3 = 0.;
	for ( unsigned int i=0; i<dim; ++i )
    	for ( unsigned int j=0; j<dim; ++j )
        	for ( unsigned int k=0; k<dim; ++k )
            	for ( unsigned int l=0; l<dim; ++l )
            	{
            		error_dfor += std::abs(D[i][j][k][l]-D_for1[i][j][k][l]);
            		error_dfor2 += std::abs(D[i][j][k][l]-D_for2[i][j][k][l]);
//            		error_for3 += std::abs(Ma_x_Ma[i][j][k][l]-Ma_x_Ma_for3[i][j][k][l]);
            	}
    
	std::cout << std::endl;
    std::cout << "D=      " << D << std::endl;
    std::cout << "D_for1= " << D_for1 << " with error=" << error_dfor << std::endl;
    std::cout << "D_for2= " << D_for2 << " with error=" << error_dfor2 << std::endl;
}

//----------------------------------------------------