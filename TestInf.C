/* This example with infinite elements computes a
*  prototypical problem of quantum mechanics:
*  A particle in a box with finite walls.
*  
* To see the effect of infinite elements, the calculation
* is performed without infinite elements as well and the solutions
* are compared.
*/
#include <string.h>
#include <complex.h>
// libMesh include files.
#include "libmesh/getpot.h" // for input-argument parsing
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
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

// Bring in everything from the libMesh namespace
using namespace libMesh;

//prototypes of functions needed to set-up the system:
void assemble_SchroedingerEquation(libMesh::EquationSystems & , const std::string &);
void get_dirichlet_dofs(EquationSystems &, const std::string & , std::set<unsigned int>&);
void tetrahedralise_sphere(UnstructuredMesh& mesh, std::vector<Node> geometry, std::string creator, Real r, int NrBall, Real VolConst, Real L, unsigned int N);
void  mesh_write(EquationSystems& equation_systems);
void  solution_write(EquationSystems& equation_systems, unsigned int i);

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
   //  INFINITE ELEMENTS ARE IMPLIMENTED ONLY FOR 3 DIMENSIONS AT THE MOMENT.
   libmesh_example_requires(3 <= LIBMESH_DIM, "3D support");
   
   int dim = 3;
   // Create two different meshes; to one of them infinite elements will be added
   // and finally the results of the calculations can be compared.
   Mesh mesh(init.comm(), dim);
   Mesh inf_mesh(init.comm(), dim);
   
   // overwrite the meshes with a spherical grid with radius 1.7
   Real E = cl("Energy", 0.0);
   Real r=cl("radius", 20.);
   //int NrBall=cl("points", 50);
   //Real VolConst= cl("maxVol", 1./(32.*sqrt(E*E*E)) );
   //Real L=cl("bending", 2.);
   //int N=cl("circles", 5);
   int maxiter=cl("maxiter", 700);
   std::vector<Node> geometry(1);
   geometry[0]=Point(0,0,0);
   std::string mesh_geom = cl("mesh_geom", "sphere");
   //tetrahedralise_sphere(mesh, geometry, mesh_geom, r, NrBall, VolConst, L, N);
 //   MeshTools::Generation::build_cube (mesh, 7, 7, 7,
 //                                         -r, r, -r, r,
 //                                         -r, r, HEX8);
 //   MeshTools::Generation::build_cube (inf_mesh, 7, 7, 7,
 //                                         -r, r, -r, r,
 //                                         -r, r, HEX8);
    MeshTools::Generation::build_sphere (mesh, r, 3, HEX8,
                                          2, true);
    MeshTools::Generation::build_sphere (inf_mesh, r, 3, HEX8,
                                          2, true);


   //In case of infinite elements, they are added now by respective interface
   InfElemBuilder builder(inf_mesh);
   builder.build_inf_elem(true);
   
   // Reassign subdomain_id() of all infinite elements.
   // Otherwise, the exodus-api will fail.
   MeshBase::element_iterator       elem_it  = inf_mesh.elements_begin();
   const MeshBase::element_iterator elem_end = inf_mesh.elements_end();
   for (; elem_it != elem_end; ++elem_it){
      Elem* elem = *elem_it;
      if(elem->infinite()){
         elem->subdomain_id() = 1;
      }
   }

   // Create an equation systems object for both systems.
   // Infinite elements don't need Condensed system in principle.
   EquationSystems finite_eq_sys (mesh);
   EquationSystems infinite_eq_sys (inf_mesh);
   
   CondensedEigenSystem & finite_eig_sys = finite_eq_sys.add_system<CondensedEigenSystem> ("EigenSE");
   CondensedEigenSystem & infinite_eig_sys = infinite_eq_sys.add_system<CondensedEigenSystem> ("EigenSE");

   // will be approximated using first-order approximation.
   finite_eig_sys.add_variable("phi", FIRST);
   infinite_eig_sys.add_variable("phi2", FIRST);

   // assign the ground state energy. For infinite elements, this guess should be 
   // good; otherwise the long-range limit will be wrong.
   finite_eq_sys.parameters.set<Number>("gsE")=E;
   infinite_eq_sys.parameters.set<Number>("gsE")=E;
   
   // set numerical parameters for SLEPC on how to solve the system.
   finite_eig_sys.eigen_solver->set_eigensolver_type(KRYLOVSCHUR); // this is default
   infinite_eig_sys.eigen_solver->set_eigensolver_type(KRYLOVSCHUR);
   //finite_eig_sys.eigen_solver->set_position_of_spectrum(SMALLEST_REAL);
   //infinite_eig_sys.eigen_solver->set_position_of_spectrum(SMALLEST_REAL);
   finite_eig_sys.eigen_solver->set_position_of_spectrum(SMALLEST_MAGNITUDE);
   infinite_eig_sys.eigen_solver->set_position_of_spectrum(SMALLEST_MAGNITUDE);

   //set number of eigen values ( \p nev) and number of 
   // basis vectors \p ncv for the solution.
   //Note that ncv >= nev must hold and ncv >= 2*nev is recommended.
   finite_eq_sys.parameters.set<unsigned int>("eigenpairs")    = nev;
   finite_eq_sys.parameters.set<unsigned int>("basis vectors") = nev*3+4;
   infinite_eq_sys.parameters.set<unsigned int>("eigenpairs")    = nev;
   infinite_eq_sys.parameters.set<unsigned int>("basis vectors") = nev*3+4;
 
   // attach the name of the function that assembles the matrix equation:
   finite_eig_sys.attach_assemble_function (assemble_SchroedingerEquation);
   infinite_eig_sys.attach_assemble_function (assemble_SchroedingerEquation);
  
   // important to set the system to be generalised hermitian eigen problem.
   // By default it is HEP and so _matrix_B is not available.
   finite_eig_sys.set_eigenproblem_type(GHEP);
   infinite_eig_sys.set_eigenproblem_type(GHEP);
   
   // Set the solver tolerance and the maximum number of iterations.
   finite_eq_sys.parameters.set<Real> ("linear solver tolerance") = pow(TOLERANCE, 5./3.);
   finite_eq_sys.parameters.set<unsigned int>("linear solver maximum iterations") = maxiter;
   infinite_eq_sys.parameters.set<Real> ("linear solver tolerance") = pow(TOLERANCE, 5./3.);
   infinite_eq_sys.parameters.set<unsigned int>("linear solver maximum iterations") = maxiter;

   
   // Initialize the data structures for the equation system.
   finite_eq_sys.init();
   infinite_eq_sys.init();
   
   // if no infinite elements are used, we need some boundary condition.
   // Here: force the wave function to vanish at all boarders.
   std::set<unsigned int> dirichlet_dof_ids;
   get_dirichlet_dofs(finite_eq_sys, "EigenSE" ,dirichlet_dof_ids);
   finite_eig_sys.initialize_condensed_dofs(dirichlet_dof_ids);

   // Solve both systems. In this function, the assemble-functions are called.
   finite_eig_sys.solve();
   infinite_eig_sys.solve();
   
   // get number of converged eigenpairs
   unsigned int fnconv = finite_eig_sys.get_n_converged();
   unsigned int inconv = infinite_eig_sys.get_n_converged();

   out << "Number of converged eigenpairs: " << fnconv << "\n";
   out << "Number of converged eigenpairs: " << inconv << "\n" << std::endl;

   #ifdef LIBMESH_HAVE_EXODUS_API
       // Write the eigen vector to file and the eigenvalues to libMesh::out.
       out<<"energy of state   (without infinite)      (with infinite)"<<std::endl;
       for(unsigned int i=0; i<std::min(inconv,fnconv); i++){
          std::pair<Real,Real> eigpair = finite_eig_sys.get_eigenpair(i);
          out<<"                   ";
          out<<eigpair.first+infinite_eq_sys.parameters.set<Number>("gsE");
          std::ostringstream eigenvector_output_name;
          eigenvector_output_name<< i <<".e" ;
          ExodusII_IO (mesh).write_equation_systems ( eigenvector_output_name.str(), finite_eq_sys);

          std::ostringstream ieigenvector_output_name;
          eigpair = infinite_eig_sys.get_eigenpair(i);
          out<<"        "<<eigpair.first+finite_eq_sys.parameters.set<Number>("gsE")<<std::endl;
          ieigenvector_output_name<< i <<"_inf.e" ;
          ExodusII_IO (inf_mesh).write_equation_systems (ieigenvector_output_name.str(), infinite_eq_sys);
       }
   #endif // #ifdef LIBMESH_HAVE_EXODUS_API
   //infinite_eig_sys.get_eigenpair(0);
   //mesh_write(finite_eq_sys);
   //mesh_write(infinite_eq_sys);
   for(unsigned int i=0; i<std::min(inconv,fnconv); i++){
      out<<std::endl<<std::endl;
      solution_write(infinite_eq_sys, i);
      out<<std::endl<<std::endl;
   }

   // All done.
   return 0;
#endif // LIBMESH_HAVE_SLEPC
}

void get_dirichlet_dofs(EquationSystems & es, const std::string & system_name, std::set<unsigned int>& dirichlet_dof_ids){
   // This function is just copied from example...
   dirichlet_dof_ids.clear();
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
   const DofMap& dof_map = eigen_system.get_dof_map();
   
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
   MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
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
      {
        // All boundary dofs are Dirichlet dofs in this case
        for (unsigned int s=0; s< elem->n_sides(); s++){
            if (elem->neighbor(s) == NULL){
               std::vector<unsigned int> side_dofs;
               FEInterface::dofs_on_side(elem, dim, fe_type, s, side_dofs);

               for(unsigned int ii=0; ii<side_dofs.size(); ii++){
                  dirichlet_dof_ids.insert(dof_indices[side_dofs[ii]]);
               }
            }
         }
      }
   } // end of element loop
   
   /**
   * All done!
   */
   return;
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
   //AutoPtr<FEBase> inf_fe (FEBase::build_InfFE(dim,fe_type));
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
   
   // A  Gauss quadrature rule for numerical integration.
   // Use the default quadrature order.
   QGauss qrule (dim, fe_type.default_quadrature_order());
      
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);
      
   libMesh::Number co0_5= 0.5;
   libMesh::Number co2= 2.;
   Number E=es.parameters.get<Number>("gsE");  
   //libMesh::Number k=omega; //divided by c which is 0 in atomic units.
   // -->ik = -i*k => for neg. energy: exp(-i*sqrt(2E)*mu(x))= exp(-sqrt(2|E|)*mu(x)) ==> expon. decay in function.
   libMesh::Number ik=sqrt(co2*E)*(std::complex<double>)_Complex_I; // -->try this for now...
   if (real(E)<0){ // E<0:
     ik=sqrt(-co2*E); // this gives exponential decay .
   }
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
            potval=-2.0;
         else
            potval=0.0; // is assumed to be close enough to infinity
         for (unsigned int i=0; i<n_sf; i++){
            for (unsigned int j=0; j<n_sf; j++){
               // this is changed here due the Petrov-Galerkin scheme. and works with finite and infinite elements.
               Se(i,j) += JxW[qp]*weight[qp]*phi[i][qp]*phi[j][qp];
               temp= dweight[qp]*phi[i][qp]*(dphi[j][qp]-ik*dphase[qp]*phi[j][qp])+
                     weight[qp]*(dphi[j][qp]*dphi[i][qp]-ik*ik*dphase[qp]*dphase[qp]*phi[i][qp]*phi[j][qp]+
                     ik*dphase[qp]*(phi[i][qp]*dphi[j][qp]-phi[j][qp]*dphi[i][qp]));
               H(i,j) += JxW[qp]*(co0_5*temp + (potval- E)*weight[qp]*phi[i][qp]*phi[j][qp]);
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

void  mesh_write(EquationSystems& equation_systems){
   CondensedEigenSystem & system = equation_systems.get_system<CondensedEigenSystem> ("EigenSE");
   const MeshBase & mesh = equation_systems.get_mesh();
   const DofMap & dof_map = system.get_dof_map();

   UniquePtr<NumericVector<Number> > solution_vect = 
        NumericVector<Number>::build(equation_systems.comm());
   solution_vect->init((*system.solution).size(), true, SERIAL);
   (*system.solution).localize(* solution_vect);
   
   const FEType & fe_type = dof_map.variable_type(0);
   UniquePtr<FEBase> fe (FEBase::build(3, fe_type));
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(3, fe_type));
   FEBase * cfe = libmesh_nullptr;
   QGauss qrule (3, SECOND);
   std::vector<dof_id_type> dof_indices;
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);
   
   MeshBase::const_element_iterator           el = mesh.active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
   for ( ; el != end_el; ++el){
      const Elem * elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      if (elem->infinite())
         cfe = inf_fe.get();
      else
          cfe = fe.get();
      const std::vector<Point>& q_point = cfe->get_xyz();
      cfe->reinit(elem);
      unsigned int max_qp = cfe->n_quadrature_points();
      for (unsigned int qp=0; qp<max_qp; qp++){
         //out<<elem->infinite()<<" ";

         //print q_point;
         out<<q_point[qp](0)<<"  ";
         out<<q_point[qp](1)<<"  ";
         out<<q_point[qp](2)<<"     ";
         Number soln=0;

         Point map_point=FEInterface::inverse_map(3, fe_type, elem, q_point[qp], TOLERANCE, true); 
         FEComputeData data(equation_systems, map_point); 
         FEInterface::compute_data(3, fe_type, elem, data);
         const unsigned int n_sf = cfe->n_shape_functions();

         //print solution value at that point.
         for (unsigned int i=0; i<n_sf; i++){
            soln+=(*solution_vect)(dof_indices[i])*data.shape[i]; // hoping the order is same in shape and dof_indices.
            //out<<std::endl<<"    "<<(*solution_vect)(dof_indices[i]);
            //out<<"    "<<data.shape[i]<<std::endl;
         }
         out<<std::real(soln)<<"  "<<std::imag(soln)<<std::endl;
      }
   }
}

void  solution_write(EquationSystems& equation_systems, unsigned int i){
   CondensedEigenSystem & system = equation_systems.get_system<CondensedEigenSystem> ("EigenSE");
   const MeshBase & mesh = equation_systems.get_mesh();
   const DofMap & dof_map = system.get_dof_map();

   // copy the i-th solution vector to system.solution.
   system.get_eigenpair(i);

   UniquePtr<NumericVector<Number> > solution_vect = 
        NumericVector<Number>::build(equation_systems.comm());

   solution_vect->init((*system.solution).size(), true, SERIAL);
   (*system.solution).localize(* solution_vect);
   
   const FEType & fe_type = dof_map.variable_type(0);
   UniquePtr<FEBase> fe (FEBase::build(3, fe_type));
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(3, fe_type));
   FEBase * cfe = libmesh_nullptr;
   QGauss qrule (3, SECOND);
   std::vector<dof_id_type> dof_indices;
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);
   
   MeshBase::const_element_iterator           el = mesh.active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
   for ( ; el != end_el; ++el){
      const Elem * elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      if (elem->infinite())
          cfe = inf_fe.get();
      else
          cfe = fe.get();
      const std::vector<Point>& q_point = cfe->get_xyz();
      cfe->reinit(elem);
      unsigned int max_qp = cfe->n_quadrature_points();
      for (unsigned int qp=0; qp<max_qp; qp++){
         //out<<elem->infinite()<<" ";

         //print q_point;
         out<<q_point[qp](0)<<"  ";
         out<<q_point[qp](1)<<"  ";
         out<<q_point[qp](2)<<"     ";
         Number soln=0;

         Point map_point=FEInterface::inverse_map(3, fe_type, elem, q_point[qp], TOLERANCE, true); 
         FEComputeData data(equation_systems, map_point); 
         FEInterface::compute_data(3, fe_type, elem, data);
         const unsigned int n_sf = cfe->n_shape_functions();

         //print solution value at that point.
         for (unsigned int i=0; i<n_sf; i++){
            soln+=(*solution_vect)(dof_indices[i])*data.shape[i]; // hoping the order is same in shape and dof_indices.
            //out<<std::endl<<"    "<<(*solution_vect)(dof_indices[i]);
            //out<<"    "<<data.shape[i]<<std::endl;
         }
         out<<std::real(soln)<<"  "<<std::imag(soln)<<std::endl;
      }
   }
}
