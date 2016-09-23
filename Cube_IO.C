#include <complex.h>
#include <iostream>
#include <fstream>
// libMesh include files.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/tree.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/point_locator_tree.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/elem.h"
#include "libmesh/eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/fe.h"
#include "libmesh/dof_map.h"
#include "libmesh/inf_fe.h"
#include "libmesh/fe_interface.h" 
#include "libmesh/fe_compute_data.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// write the solutions values at points in a cube.
void cube_io(EquationSystems& es, std::vector<Node> geom, std::string output, std::string SysName){
   //CondensedEigenSystem & system = es.get_system<CondensedEigenSystem> ("EigenSE"); // --> how to generalise??
   System & system = es.get_system<System> (SysName); 
   const MeshBase & mesh = es.get_mesh();
   const DofMap & dof_map = system.get_dof_map();
   
   UniquePtr<NumericVector<Number> > solution_vect = 
        NumericVector<Number>::build(es.comm());

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
   re_out<<SysName<<std::endl<<std::endl; // print first two lines: comments
   im_out<<SysName<<std::endl<<std::endl; 
   abs_out<<SysName<<std::endl<<std::endl;

   re_out<<std::setw(5)<<"  "<<geom.size();
   im_out<<std::setw(5)<<"  "<<geom.size();
   abs_out<<std::setw(5)<<"  "<<geom.size();
   // where do I start?
   Point mol_center;
   Point min(0,0,0);
   Point max(0,0,0);
   for(unsigned int i=0; i<geom.size(); i++){
      mol_center(0)+=geom[i](0);
      mol_center(1)+=geom[i](1);
      mol_center(2)+=geom[i](2);
      if(geom[i](0)<min(0))
         min(0)=geom[i](0);
      if(geom[i](0)>max(0))
         max(0)=geom[i](0);
      if(geom[i](1)<min(1))
         min(1)=geom[i](1);
      if(geom[i](1)>max(1))
         max(2)=geom[i](1);
      if(geom[i](2)<min(2))
         min(2)=geom[i](2);
      if(geom[i](2)>max(2))
         max(2)=geom[i](2);
   }
   mol_center(0)=mol_center(0)/geom.size();
   mol_center(1)=mol_center(1)/geom.size();
   mol_center(2)=mol_center(2)/geom.size();

   Real dx=0.5;
   Real dy=0.5;
   Real dz=0.5;
   unsigned int nx=60+(max(0)-min(0))*2;
   unsigned int ny=60+(max(1)-min(1))*2;
   unsigned int nz=60+(max(2)-min(2))*2;

   Point start(mol_center(0)-dx*nx/2,
               mol_center(1)-dy*ny/2,
               mol_center(2)-dz*nz/2);

   re_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(0);
   re_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(1);
   re_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(2)<<std::endl;
   im_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(0);
   im_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(1);
   im_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(2)<<std::endl;
   abs_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(0);
   abs_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(1);
   abs_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(2)<<std::endl;
   // print # points per axis and step in Cartesian Coordinates:

   re_out<<std::setw(5)<<nx;
   re_out<<std::setw(12)<<std::setprecision(5)<<" \t "<<dx<<" \t 0.00000 \t 0.00000"<<std::endl;
   re_out<<std::setw(5)<<ny;
   re_out<<std::setw(12)<<std::setprecision(5)<<" \t 0.00000 \t "<<dy<<" \t 0.00000"<<std::endl;
   re_out<<std::setw(5)<<nz;
   re_out<<std::setw(12)<<std::setprecision(5)<<" \t 0.00000 \t 0.00000 \t "<<dz<<std::endl;
   im_out<<std::setw(5)<<nx;
   im_out<<std::setw(12)<<std::setprecision(5)<<" \t "<<dx<<" \t 0.00000 \t 0.00000"<<std::endl;
   im_out<<std::setw(5)<<ny;
   im_out<<std::setw(12)<<std::setprecision(5)<<" \t 0.00000 \t "<<dy<<" \t 0.00000"<<std::endl;
   im_out<<std::setw(5)<<nz;
   im_out<<std::setw(12)<<std::setprecision(5)<<" \t 0.00000 \t 0.00000 \t "<<dz<<std::endl;
   abs_out<<std::setw(5)<<nx;
   abs_out<<std::setw(12)<<std::setprecision(5)<<" \t "<<dx<<" \t 0.00000 \t 0.00000"<<std::endl;
   abs_out<<std::setw(5)<<ny;
   abs_out<<std::setw(12)<<std::setprecision(5)<<" \t 0.00000 \t "<<dy<<" \t 0.00000"<<std::endl;
   abs_out<<std::setw(5)<<nz;
   abs_out<<std::setw(12)<<std::setprecision(5)<<" \t 0.00000 \t 0.00000 \t "<<dz<<std::endl;

   for(unsigned int i=0; i<geom.size(); i++){
      re_out<<" "<<std::setw(5)<<geom[i].id()<<"\t";
      im_out<<" "<<std::setw(5)<<geom[i].id()<<"\t";
      abs_out<<" "<<std::setw(5)<<geom[i].id()<<"\t";
      //out<<std::setw(5)<<"1.00000"<<"\t";
      re_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      re_out<<" "<<std::setw(12)<<std::setprecision(6)<<geom[i](0)<<"\t";
      re_out<<" "<<std::setw(12)<<std::setprecision(6)<<geom[i](1)<<"\t";
      re_out<<" "<<std::setw(12)<<std::setprecision(6)<<geom[i](2)<<"\n";
      im_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      im_out<<" "<<std::setw(12)<<std::setprecision(6)<<geom[i](0)<<"\t";
      im_out<<" "<<std::setw(12)<<std::setprecision(6)<<geom[i](1)<<"\t";
      im_out<<" "<<std::setw(12)<<std::setprecision(6)<<geom[i](2)<<"\n";
      abs_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      abs_out<<" "<<std::setw(12)<<std::setprecision(6)<<geom[i](0)<<"\t";
      abs_out<<" "<<std::setw(12)<<std::setprecision(6)<<geom[i](1)<<"\t";
      abs_out<<" "<<std::setw(12)<<std::setprecision(6)<<geom[i](2)<<"\n";
   }

   unsigned int ix, iy, iz;
   PointLocatorTree pt_lctr(mesh);
   unsigned int num_line=0;
   for (ix=0;ix<nx;ix++) {
      for (iy=0;iy<ny;iy++) {
         for (iz=0;iz<nz;iz++) {

            num_line++;
            Point q_point(start(0)+(Real)ix*dx,
                          start(1)+(Real)iy*dy,
                          start(2)+(Real)iz*dz);
            
            const Elem * elem=pt_lctr(q_point);
            if(elem==NULL){
               abs_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<0.0;
               im_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<0.0;
               re_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<0.0;
            }
            else{

               dof_map.dof_indices (elem, dof_indices);
   
               Point map_point=FEInterface::inverse_map(3, fe_type, elem, q_point, TOLERANCE, true); 
               if (elem->infinite()){
                  out<<map_point<<std::endl;
               }
               FEComputeData data(es, map_point); 
               if (elem->infinite()){
                  out<<"point: "<<q_point<<"   ";
                  out<<map_point(2)<<std::endl;
               }
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
               re_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::real(soln);
               im_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::imag(soln);
               abs_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::abs(soln);
            }

            if (num_line == 6){
               re_out<<std::endl;
               im_out<<std::endl;
               abs_out<<std::endl;
               num_line=0;
            }
         }
         if (num_line>0){
            re_out<<std::endl;
            im_out<<std::endl;
            abs_out<<std::endl;
         }
      }
   }
}

// write the solutions values at quadrature points
void  solution_write(EquationSystems& equation_systems, std::string filename, std::string SysName){
   System & system = equation_systems.get_system<System> (SysName); 
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


   std::ofstream out(filename);
   out<<std::endl<<std::endl;
   
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
   out<<std::endl<<std::endl;
}
