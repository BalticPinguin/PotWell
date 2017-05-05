/* This example with infinite elements computes a
*  prototypical problem of quantum mechanics:
*  A particle in a box with finite walls.
*  
* To see the effect of infinite elements, the calculation
* is performed without infinite elements as well and the solutions
* are compared.
*/
#include <string.h>
#include <iostream>
#include <fstream>
#include <complex.h>
// libMesh include files.
#include "libmesh/getpot.h" // for input-argument parsing
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/fe_interface.h" // for dirichlet boundary conditions
#include "libmesh/error_vector.h" // for dirichlet boundary conditions
#include "libmesh/explicit_system.h"
// for infinite elements:
#include "libmesh/inf_fe.h"
#include "libmesh/inf_elem_builder.h"
// for dirichlet boundary conditions
#include "libmesh/fe_interface.h" 
#include "libmesh/fe_compute_data.h"
#include "libmesh/error_vector.h" 

// for finding element for point
#include "libmesh/point_locator_tree.h"

// for the SlepcSolverConfiguration
#include "libmesh/solver_configuration.h"
#include "libmesh/slepc_eigen_solver.h"


EXTERN_C_FOR_SLEPC_BEGIN
# include <slepceps.h>
EXTERN_C_FOR_SLEPC_END

/**
 * Defines an \p enum for spectral tronsformations
 * applied before solving the (generalised) eigenproblem
 */
enum SpectralTransform {SHIFT=0,
                        SINVERT,
                        CAYLEY,

                        INVALID_ST
};

class SlepcSolverConfiguration : public libMesh::SolverConfiguration
{
public:

   SlepcSolverConfiguration( libMesh::SlepcEigenSolver<libMesh::Number> & slepc_eigen_solver):
        _slepc_solver(slepc_eigen_solver),
        _st(INVALID_ST)
   {}
   
   ~SlepcSolverConfiguration() {}

   virtual void configure_solver() override;

   void SetST(SpectralTransform st)
   { _st=st;}
   
private:
   // The linear solver object that we are configuring
   libMesh::SlepcEigenSolver<libMesh::Number>& _slepc_solver;
   SpectralTransform _st;
   //ST st;

};

// Bring in everything from the libMesh namespace
using namespace libMesh;

//prototypes of functions needed to set-up the system:
void assemble_SchroedingerEquation(libMesh::EquationSystems & , const std::string &);
void line_print(EquationSystems& es, std::string output, std::string SysName);

int main (int argc, char** argv){
   // Initialize libMesh and the dependent libraries.
   LibMeshInit init (argc, argv);
   GetPot cl(argv[1]);

   // Skip SLEPc examples on a non-SLEPc libMesh build
   #ifndef LIBMESH_HAVE_SLEPC
     libmesh_example_requires(false, "--enable-slepc");
     }
   #else
   #ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
   libmesh_example_requires(false, "--enable-ifem");
   #endif
   
   #ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
    // SLEPc currently gives us a nasty crash with Real==float
    // FIXME: I am not sure, thi is still valid. That command is more than 10 years old...
    libmesh_example_requires(false, "--disable-singleprecision");
   #endif
   const unsigned int nev = cl("nev",5);
   #ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
   libmesh_example_requires(false, "--enable-ifem");
   #endif
   
   // Check for proper usage.
   if (argc < 2)
      libmesh_error_msg("\nUsage: " << argv[0] << " <input-filename>");
   // Tell the user what we are doing.
   else {
      std::cout << "Running " << argv[0];

      for (int i=1; i<argc; i++)
          std::cout << " " << argv[i];
      std::cout << std::endl << std::endl;
   }

   // Skip this example if libMesh was compiled with <3 dimensions.
   // INFINITE ELEMENTS ARE IMPLEMENTED ONLY FOR 3 DIMENSIONS AT THE MOMENT.
   libmesh_example_requires(3 <= LIBMESH_DIM, "3D support");
   
   int dim = 3;
   // Create two different meshes; to one of them infinite elements will be added
   // and finally the results of the calculations can be compared.
   Mesh mesh(init.comm(), dim);
   
   // overwrite the meshes with a spherical grid with radius 1.7
   Real E = cl("Energy", 0.1);
   Real r=cl("radius", 2.);
   unsigned int maxiter=cl("maxiter", 700);

   //MeshTools::Generation::build_sphere (mesh, r, 3, HEX27, 3, true);
   //MeshTools::Generation::build_sphere (mesh, r, 1, HEX8, 1, true);
   MeshTools::Generation::build_sphere (mesh, r, 2, HEX8, 2, true);

   //FEType fe_type(SECOND, LAGRANGE, FIFTH, JACOBI_20_00, CARTESIAN);
   FEType fe_type(FIRST, LAGRANGE, FIFTH, JACOBI_20_00, CARTESIAN);
   //In case of infinite elements, they are added now by respective interface

   InfElemBuilder builder(mesh);
   builder.build_inf_elem(true);
   
   // Reassign subdomain_id() of all infinite elements.
   // Otherwise, the exodus-api will fail.
   MeshBase::element_iterator       elem_it  = mesh.elements_begin();
   const MeshBase::element_iterator elem_end = mesh.elements_end();
   for (; elem_it != elem_end; ++elem_it){
      Elem* elem = *elem_it;
      if(elem->infinite()){
         elem->subdomain_id() = 1;
      }
   }

   // find the neighbours; for correct linking the two areas
   mesh.find_neighbors();

   // convert all elements from 1-st order to second order elements.
   //mesh.all_second_order();

   // Create an equation systems object for both systems.
   // Infinite elements don't need Condensed system in principle.
   EquationSystems eq_sys (mesh);
   
   CondensedEigenSystem & eig_sys = eq_sys.add_system<CondensedEigenSystem> ("EigenSE");

   // will be approximated using first-order approximation.
   eig_sys.add_variable("phi", fe_type);

   // assign the ground state energy. For infinite elements, this guess should be 
   // good; otherwise the long-range limit will be wrong.
   eq_sys.parameters.set<Number>("gsE")=E;
   
   // set numerical parameters for SLEPC on how to solve the system.
   eig_sys.eigen_solver->set_eigensolver_type(KRYLOVSCHUR); // this is default

   eig_sys.eigen_solver->set_position_of_spectrum(E);
   
   SlepcEigenSolver<Number>* solver = 
                 libmesh_cast_ptr<SlepcEigenSolver<Number>* >( &(*eig_sys.eigen_solver) );

   SlepcSolverConfiguration ConfigSolver( *solver);

   // set the spectral transformation:
   ConfigSolver.SetST(SINVERT);
   //ConfigSolver.SetST(CAYLEY);
   //ConfigSolver.SetST(SHIFT); // this is default
   solver ->set_solver_configuration(ConfigSolver);

   eq_sys.parameters.set<Real>("radius")    = r;
   //set number of eigen values ( \p nev) and number of 
   // basis vectors \p ncv for the solution.
   //Note that ncv >= nev must hold and ncv >= 2*nev is recommended.
   eq_sys.parameters.set<unsigned int>("eigenpairs")    = nev;
   eq_sys.parameters.set<unsigned int>("basis vectors") = nev*3+4;
 
   // attach the name of the function that assembles the matrix equation:
   eig_sys.attach_assemble_function (assemble_SchroedingerEquation);
  
   // important to set the system to be generalised nonhermitian eigen problem.
   // By default it is HEP and so _matrix_B is not available.
   eig_sys.set_eigenproblem_type(GNHEP);
   
   // Set the solver tolerance and the maximum number of iterations.
   eq_sys.parameters.set<Real> ("linear solver tolerance") = pow(TOLERANCE, 5./3.);
   eq_sys.parameters.set<unsigned int>("linear solver maximum iterations") = maxiter;

   // Initialize the data structures for the equation system.
   eq_sys.init();

   // Solve system. In this function, the assemble-functions are called.
   eig_sys.solve();
   
   // get number of converged eigenpairs
   unsigned int nconv = eig_sys.get_n_converged();

   out << "Number of converged eigenpairs: " << nconv << "\n";

   // Write the eigen vector to file and the eigenvalues to libMesh::out.
   for(unsigned int i=0; i<nconv; i++){
      std::pair<Real,Real> eigpair = eig_sys.get_eigenpair(i);
      eigpair = eig_sys.get_eigenpair(i);
      out<<"        "<<eigpair.first<<std::endl;
      eq_sys.parameters.set<Real>("current frequency")=eq_sys.parameters.get<Real>("speed")*sqrt(std::abs(eigpair.first)*2.)/(2*pi);
      
      std::ostringstream file;
      file<<"infini_"<<i<<".txt";
      // print the solution along the x-coordinate
      line_print(eq_sys, file.str(), "EigenSE");
   }

   // All done.
   return 0;
#endif // LIBMESH_HAVE_SLEPC
}

void assemble_SchroedingerEquation(libMesh::EquationSystems & es, const std::string & system_name){
   // It is a good idea to make sure we are assembling
   // the proper system.
   libmesh_assert_equal_to (system_name, "EigenSE");
   // Get a constant reference to the mesh object.
   const MeshBase& mesh = es.get_mesh();
   // The dimension that we are running.
   const unsigned int dim = mesh.mesh_dimension();
      
   // Get a reference to our system.
   CondensedEigenSystem & eigen_system = es.get_system<CondensedEigenSystem> (system_name);

   // Get a constant reference to the Finite Element type
   // for the first (and only) variable in the system.
   FEType fe_type = eigen_system.get_dof_map().variable_type(0);
      
   // A reference to the system matrix
   SparseMatrix<Number>&  matrix_A = *eigen_system.matrix_A;
   SparseMatrix<Number>&  matrix_B = *eigen_system.matrix_B;
   // Build a Finite Element object of the specified type.  Since the
   // \p FEBase::build() member dynamically creates memory we will
   // store the object as an \p UniquePtr<FEBase>.  This can be thought
   // of as a pointer that will clean up after itself.
   UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));  // here, try AutoPtr instead...
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
   
   // A  Gauss quadrature rule for numerical integration.
   // Use the default quadrature order.
   QGauss qrule (dim, fe_type.default_quadrature_order());
      
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);
      
   libMesh::Number co2= 2.;
   Number E=es.parameters.get<Number>("gsE");  
   //libMesh::Number k=omega; //divided by c which is 0 in atomic units.
   // -->ik = -i*k => for neg. energy: exp(-i*sqrt(2E)*mu(x))= exp(-sqrt(2|E|)*mu(x)) ==> expon. decay in function.
   libMesh::Number ik=sqrt(-co2*E); // -->try this for now...
   // set parameters for infinite elements:
   es.parameters.set<Real>("speed")=137.0359991;
   // --> it would be better if 'current frequency' could be <Number>, not <Real>.
   es.parameters.set<Real>("current frequency")=es.parameters.get<Real>("speed")*sqrt(std::abs(E)*2.)/(2*pi);

   Number potval;  
   libMesh::Number temp; // -->try this for now...
      
   // A reference to the \p DofMap object for this system.  The \p DofMap
   // object handles the index translation from node and element numbers
   // to degree of freedom numbers.
   const DofMap& dof_map = eigen_system.get_dof_map();
      
   // The element mass matrix and Hamiltonian
   DenseMatrix<Number> Se;
   DenseMatrix<Number> H;
   
   // This vector will hold the degree of freedom indices for
   // the element.  These define where in the global system
   // the element degrees of freedom get mapped.
   std::vector<dof_id_type> dof_indices;
      
   // Now we will loop over all the elements in the mesh that
   // live on the local processor. We will compute the element
   // matrix and right-hand-side contribution.  In case users
   // later modify this program to include refinement, we will
   // be safe and will only consider the active elements;
   // hence we use a variant of the \p active_elem_iterator.
   MeshBase::const_element_iterator       el  = mesh.active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
      
   for ( ; el != end_el; ++el){
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // unifyging finite and infinite elements
      FEBase * cfe = libmesh_nullptr;

      if (elem->infinite()){
         // We have an infinite element.  Let \p cfe point
         // to our \p InfFE object.  This is handled through
         // an UniquePtr.  Through the \p UniquePtr::get() we "borrow"
         // the pointer, while the \p  UniquePtr \p inf_fe is
         // still in charge of memory management.
         cfe = inf_fe.get();
      }
      else{
        cfe = fe.get();
      }
     
      // The element Jacobian * quadrature weight at each integration point.
      const std::vector<Real>& JxW = cfe->get_JxW();

      // The element shape functions evaluated at the quadrature points.
      const std::vector<std::vector<Real> >& phi = cfe->get_phi();
      const std::vector<std::vector<RealGradient> >& dphi = cfe->get_dphi();
      const std::vector<Point>& q_point = cfe->get_xyz();
      // get extra data needed for infinite elements
      const std::vector<RealGradient>& dphase = cfe->get_dphase();
      const std::vector<Real>& weight = cfe->get_Sobolev_weight(); // in publication called D
      const std::vector<RealGradient>& dweight = cfe->get_Sobolev_dweight();
   
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      cfe->reinit (elem);
   
      // Zero the element matrices and rhs before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Se.resize (dof_indices.size(), dof_indices.size());
      H.resize (dof_indices.size(), dof_indices.size());

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      //For infinite elements, the number of quadrature points is asked and than looped over; works for finite elements as well.
      unsigned int max_qp = cfe->n_quadrature_points();
      for (unsigned int qp=0; qp<max_qp; qp++){
         // Now, get number of shape functions that are nonzero at this point::
         unsigned int n_sf = cfe->n_shape_functions();
         // loop over them:
         if (q_point[qp].size()<1.)
            potval=-0.5;
         else
            potval=0.0; // is assumed to be close enough to infinity
         for (unsigned int i=0; i<n_sf; i++){
            for (unsigned int j=0; j<n_sf; j++){
               // this is changed here due the Petrov-Galerkin scheme. and works with finite and infinite elements.
               Se(i,j) += JxW[qp]*weight[qp]*phi[i][qp]*phi[j][qp];
               temp= dweight[qp]*phi[i][qp]*(dphi[j][qp]-ik*dphase[qp]*phi[j][qp])+
                     weight[qp]*(dphi[j][qp]*dphi[i][qp]-ik*ik*dphase[qp]*dphase[qp]*phi[i][qp]*phi[j][qp]+

                     ik*dphase[qp]*(phi[i][qp]*dphi[j][qp]-phi[j][qp]*dphi[i][qp]));
               H(i,j) += JxW[qp]*(0.5*temp + potval*weight[qp]*phi[i][qp]*phi[j][qp]);
            }
         }
      }
      // On an unrefined mesh, constrain_element_matrix does
      // nothing.  If this assembly function is ever repurposed to
      // run on a refined mesh, getting the hanging node constraints
      // right will be important.  Note that, even with
      // asymmetric_constraint_rows = false, the constrained dof
      // diagonals still exist in the matrix, with diagonal entries
      // that are there to ensure non-singular matrices for linear
      // solves but which would generate positive non-physical
      // eigenvalues for eigensolves.
      dof_map.constrain_element_matrix(Se, dof_indices, false);
      dof_map.constrain_element_matrix(H, dof_indices, false);

      // Finally, simply add the element contribution to the
      // overall matrix.
      matrix_A.add_matrix (H, dof_indices);
      matrix_B.add_matrix (Se, dof_indices);

   } // end of element loop
         
   /**
   * All done!
   */
   return;
}

void SlepcSolverConfiguration::configure_solver()
{
   PetscErrorCode ierr = 0;

   // if a spectral transformation was requested
   if (_st!=INVALID_ST){
    
      // initialise the st with the default values
      //(than, change only the spectral transformation value).
      ST st;
      ierr = EPSGetST(_slepc_solver.eps(), &st);
      libmesh_assert(ierr == 0);
      //STCreate(_slepc_solver.comm().get(), &st);

      // Set it to the desired type of spectral transformation.
      // The value of the respective shift is chosen to be the target
      // specified via \p set_position_of_spectrum().
      switch (_st)
         {
         case SHIFT:
            ierr = STSetType(st, STSHIFT);
            break;
         case SINVERT:
            ierr = STSetType(st, STSINVERT);
            break;
         case CAYLEY:
      #if SLEPC_VERSION_LESS_THAN(2,2,1)
            libmesh_error_msg("SLEPc 2.2.1 is required to call CAYLEY transform.");
            break;
      #else
            ierr = STSetType(st, STCAYLEY);
            break;
      #endif
         default:
            // print a warning but do nothing more.
            break;
         }  //tell the \p EPS object which \p ST to use
      // this is not needed because it is called in the
      // in the \p EPSSetUP() anyway.
      //ierr = EPSSetST(_slepc_solver.eps(), st);

      libmesh_assert(ierr == 0);
   }
}

void line_print(EquationSystems& es, std::string output, std::string SysName){
   //CondensedEigenSystem & system = es.get_system<CondensedEigenSystem> ("EigenSE"); // --> how to generalise??
   System & system = es.get_system<System> (SysName); 
   const MeshBase & mesh = es.get_mesh();
   const DofMap & dof_map = system.get_dof_map();
   
   UniquePtr<NumericVector<Number> > solution_vect = 
        NumericVector<Number>::build(es.comm());

   solution_vect->init((*system.solution).size(), true, SERIAL);
   (*system.solution).localize(* solution_vect);
   Real r = es.parameters.get<Real>("radius");
   
   const FEType & fe_type = dof_map.variable_type(0);
   UniquePtr<FEBase> fe (FEBase::build(3, fe_type));
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(3, fe_type));
   FEBase * cfe = libmesh_nullptr;
   QGauss qrule (3, SECOND);
   std::vector<dof_id_type> dof_indices;
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);

   // set output to filename
   std::ostringstream re_output;
   re_output<<"re_"<<output;
   std::ostringstream im_output;
   im_output<<"im_"<<output;
   std::ostringstream abs_output;
   abs_output<<"abs_"<<output;
   std::ofstream im_out(re_output.str());
   std::ofstream re_out(im_output.str());
   std::ofstream abs_out(abs_output.str());

   PointLocatorTree pt_lctr(mesh);
   unsigned int num_line=0;
   Real N = 100.;
   Point q_point;
   const Real start=-2*r;
   for (int pts=1;pts<=2*N;pts++) {
      // go from -2*r to 2*r.
      q_point = Point( start+ pts*r/N*2, 0., 0.);
      num_line++;
      
      const Elem * elem=pt_lctr(q_point);
      if(elem==NULL){
         re_out<<" "<<std::setw(12)<<q_point(0);
         im_out<<" "<<std::setw(12)<<q_point(0);
         abs_out<<" "<<std::setw(12)<<q_point(0);
         abs_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<0.0<<std::endl;
         im_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<0.0<<std::endl;
         re_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<0.0<<std::endl;
      }
      else{

         dof_map.dof_indices (elem, dof_indices);
  
         Point map_point=FEInterface::inverse_map(3, fe_type, elem, q_point, TOLERANCE, true); 
         FEComputeData data(es, map_point); 
         FEInterface::compute_data(3, fe_type, elem, data);
      
         //compute solution value at that point.
         Number soln=0;
         if (elem->infinite())
            cfe = inf_fe.get();
         else
            cfe = fe.get();
         cfe->reinit(elem);
         unsigned int n_sf= cfe->n_shape_functions();
         for (unsigned int i=0; i<n_sf; i++){
            soln+=(*solution_vect)(dof_indices[i])*data.shape[i];
         }
         re_out<<" "<<std::setw(12)<<q_point(0);
         im_out<<" "<<std::setw(12)<<q_point(0);
         abs_out<<" "<<std::setw(12)<<q_point(0);
         re_out<<"  "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::real(soln)<<std::endl;
         im_out<<"  "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::imag(soln)<<std::endl;
         abs_out<<"  "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::abs(soln)<<std::endl;

      }
   }
}
