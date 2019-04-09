#include "FEM3.h"

// Class constructor for a vector field
template <int dim>
FEM<dim>::FEM (unsigned int prob)
:
fe (FE_Q<dim>(order), dim),
  dof_handler (triangulation),
  quadrature_formula(quadRule),
  face_quadrature_formula(quadRule)
{       
        
  problem = prob; //Should be 1 or 2
        
  //Nodal Solution names - this is for writing the output file
  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u");
    nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }
}

//Class destructor
template <int dim>
FEM<dim>::~FEM (){dof_handler.clear ();}


/**************************************************************************************/
/*  function C returns the components of the 4th order elasticity tensor.             */
/*        C(i,j,k,l) = lambda*1(i,j)*1(k,l) + mu*[1(i,k)*1(j,l)+1(i,l)*1(j,k)]        */
/*    where 1(x, y) is the Kronecker delta                                            */
/**************************************************************************************/

template <int dim>
double FEM<dim>::C(unsigned int i,unsigned int j,unsigned int k,unsigned int l){

  //Define the material parameters of Young's modulus and Poisson's ratio
  double E= 2.0e11,
    nu= 0.3;
  double lambda=(E*nu)/((1.+nu)*(1.-2.*nu)),
    mu=E/(2.*(1.+nu));

  return lambda*(i==j)*(k==l) + mu*((i==k)*(j==l) + (i==l)*(j==k));

}


/**************************************************************************************/
/*  Define the problem domain and generate the mesh                                   */
/**************************************************************************************/

template <int dim>
void FEM<dim>::generate_mesh(std::vector<unsigned int> numberOfElements){

  //Define the limits of your domain
    x_min = 0.,
    x_max = 1.,
    y_min = 0.,
    y_max = 1.,
    z_min = 0.,
    z_max = 1.;

  Point<dim,double> min(x_min,y_min,z_min),
    max(x_max,y_max,z_max);
  GridGenerator::subdivided_hyper_rectangle (triangulation, numberOfElements, min, max);
}


/**************************************************************************************/
/*  This function specifies and stores Dirichilet boundary conditions.                */
/*  After assembling Kglobal and Flocal, theses b.c.s will be automatically applied   */
/*  to resive matrices of the system with techniques supported by deal.II.            */
/**************************************************************************************/

template <int dim>
void FEM<dim>::define_boundary_conds(){

  const unsigned int totalDOFs = dof_handler.n_dofs(); //Total number of degrees of freedom
    unsigned int nodalDOF;
    
    //loop over global degrees of freedom 
    for(unsigned int global_dof=0; global_dof<totalDOFs; global_dof++){

      if(problem == 1) {
        // Since in both problems, u1 = u2 = u3 = 0 on specified boundary, we do not
        // need to add an additional judging condition of nodalDOF = globalDOF % dim
        if(dofLocation[global_dof][2] == z_min) { // z = 1.
          boundary_values[global_dof] = 0.;
        }
      } // if(problem == 1)

      else if(problem == 2) {
        if(dofLocation[global_dof][0] == x_min) { // x = 0.
          boundary_values[global_dof] = 0.;
        }
        else if(dofLocation[global_dof][0] == x_max) { // x = 1.
          double x2 = dofLocation[global_dof][1], x3 = dofLocation[global_dof][2];
          if(global_dof % dim == 1)
            boundary_values[global_dof] = 0.5 + (x2 - 0.5) * cos(PI / 30.) - ( x3 - 0.5) * sin(PI / 30.) - x2;
          if(global_dof % dim == 2)
            boundary_values[global_dof] = 0.5 + (x2 - 0.5) * sin(PI / 30.) + ( x3 - 0.5) * cos(PI / 30.) - x3;
        }
      } // if(problem == 2)
    }
}


/**************************************************************************************/
/*  This function sets up data structures for Finite Element Solving System.          */
/**************************************************************************************/

template <int dim>
void FEM<dim>::setup_system(){

  //Let deal.II organize degrees of freedom
  dof_handler.distribute_dofs (fe);

  //Get a vector of global degree-of-freedom x-coordinates
  MappingQ1<dim,dim> mapping;
  std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
  dofLocation.reinit(dof_handler.n_dofs(),dim);
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++){
    for(unsigned int j=0; j<dim; j++){
      dofLocation[i][j] = dof_coords[i][j];
    }
  }

  //Specify boundary condtions (call the function)
  define_boundary_conds();

  //Define the size of the global matrices and vectors
  sparsity_pattern.reinit (dof_handler.n_dofs(), 
                           dof_handler.n_dofs(),
                           dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
  K.reinit (sparsity_pattern);
  D.reinit(dof_handler.n_dofs());
  F.reinit(dof_handler.n_dofs());

  //Just some notes...
  std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;   
}


/**************************************************************************************/
/*  This function assembles K_local, F_local, K_glocal and F_global matrices/vectors. */
/*  It generally applies a quadrature approximation to estimate an itegral over the   */
/*  reference domain, rather than solve it analytically.                              */
/**************************************************************************************/

template <int dim>
void FEM<dim>::assemble_system(){

  //For volume integration/quadrature points
  FEValues<dim> fe_values (fe,
                             quadrature_formula, 
                             update_values | 
                             update_gradients | 
                             update_JxW_values);

  //For surface integration/quadrature points
  FEFaceValues<dim> fe_face_values (fe,
                                      face_quadrature_formula, 
                                      update_values | 
                                      update_quadrature_points | 
                                      update_JxW_values);

  K=0; F=0;

  const unsigned int dofs_per_elem = fe.dofs_per_cell;                      //dofs per element
  const unsigned int nodes_per_elem = GeometryInfo<dim>::vertices_per_cell;
  const unsigned int num_quad_pts = quadrature_formula.size();              //Total number of quad points in the element
  const unsigned int num_face_quad_pts = face_quadrature_formula.size();    //Total number of quad points in the face
  const unsigned int faces_per_elem = GeometryInfo<dim>::faces_per_cell;
  FullMatrix<double> Klocal (dofs_per_elem, dofs_per_elem);
  Vector<double>     Flocal (dofs_per_elem);

  std::vector<unsigned int> local_dof_indices (dofs_per_elem);              //relates local dof numbering to global dof numbering

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    fe_values.reinit(elem);

    elem->get_dof_indices (local_dof_indices);
                
    Klocal = 0.;

    //evaluate elemental stiffness matrix
    //  K^{AB}_{ik} = \integral N^A_{,j}*C_{ijkl}*N^B_{,l} dV 
    for(unsigned int q = 0; q < num_quad_pts; q++) {
      for(unsigned int A = 0; A < dofs_per_elem / dim; A++) {
        for(unsigned int i = 0; i < dim; i++) {
          for(unsigned int B = 0; B < dofs_per_elem / dim; B++) {
            for(unsigned int k = 0; k < dim; k++) {
              for(unsigned int j = 0; j < dim; j++) {
                for(unsigned int l = 0; l < dim; l++) {
                  Klocal[dim*A+i][dim*B+k] += fe_values.shape_grad(dim*A+i,q)[j] * C(i,j,k,l)
                                    * fe_values.shape_grad(dim*B+k,q)[l] * fe_values.JxW(q);
                }
              }
            }
          }
        }
      }
    } // for(q)

    Flocal = 0.;

    //Add Neumann boundary conditions here in Flocal by integrating over the appropriate surface
      Vector<double> h(dim); h=0.;
      for (unsigned int f=0; f < faces_per_elem; f++){
        //Update fe_face_values from current element and face
        fe_face_values.reinit (elem, f);

        if(elem->face(f)->center()[2] == 1){ // x3 == 1
          //To integrate over this face, loop over all face quadrature points with this single loop
          for (unsigned int q=0; q<num_face_quad_pts; ++q){
            if(problem == 1) {
              double x; //x-coordinate at the current surface quad. point
              x = fe_face_values.quadrature_point(q)[0];
              h[2] = 1.0e9 * x;
            }
            for (unsigned int A=0; A<dofs_per_elem / dim; A++){ //loop over all element nodes
              for(unsigned int i=0; i<dim; i++){ //loop over nodal dofs
                Flocal[dim*A+i] += fe_face_values.shape_value(dim*A+i, q) * h[i] * fe_face_values.JxW(q);
              }
            }
          }
        }
      }

    //Assemble local K and F into global K and F
    for(unsigned int i=0; i<dofs_per_elem; i++){
      F[local_dof_indices[i]] += Flocal[i]; 
      for(unsigned int j=0; j<dofs_per_elem; j++){
        K.add(local_dof_indices[i], local_dof_indices[j], Klocal[i][j]);
      }
    }
  } // for(elem)

  //apply Dirichlet conditions WITHOUT modifying the size of K and F global
  MatrixTools::apply_boundary_values (boundary_values, K, D, F, false);
}


/**************************************************************************************/
/*  This function solves for the torque about the x1 axis on the face x1=1            */
/*  The torque is calculated as \integral (T * arm) on a face                         */
/**************************************************************************************/

template <int dim>
double FEM<dim>::torque(){

  FEFaceValues<dim> fe_face_values (fe,
                                    face_quadrature_formula, 
                                    update_values |  
                                    update_gradients | 
                                    update_quadrature_points | 
                                    update_JxW_values);

  //Integrate the two shear stresses (multiplied by the appropriate moment arms)
  double torque = 0.;

  unsigned int faces_per_elem = GeometryInfo<dim>::faces_per_cell;
  const unsigned int                              nodes_per_elem = GeometryInfo<dim>::vertices_per_cell;
  const unsigned int dofs_per_elem = fe.dofs_per_cell; //This gives you dofs per element
  unsigned int num_face_quad_pts = face_quadrature_formula.size(); //Total number of qua
  std::vector<unsigned int> local_dof_indices (dofs_per_elem); //This relates local dof numbering to global dof numbering
  double sigma_12, sigma_13;
  FullMatrix<double> du_dx(dim,dim);
  unsigned int elemDOF;

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (), endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    //Retrieve the effective "connectivity matrix" for this element
    elem->get_dof_indices (local_dof_indices);

    //Loop over faces; integrate over the face x1=1
    for (unsigned int f=0; f < faces_per_elem; f++){
      if(elem->face(f)->center()[0] == 1){
        fe_face_values.reinit (elem, f); //Retrieve values from current element and face
        for (unsigned int q=0; q<num_face_quad_pts; ++q){
          // define the moment arms along y and z-axis
          double yarm, zarm;
          yarm = fe_face_values.quadrature_point(q)[1];
          zarm = zarm = fe_face_values.quadrature_point(q)[2];
          sigma_12 = 0.; sigma_13 = 0.;
          du_dx = 0;
          for (unsigned int A=0; A<nodes_per_elem; A++){ //loop over all element nodes
            for(unsigned int k=0; k<dim; k++){ //loop over nodal dofs
              elemDOF = dim*A+k;
              for(unsigned int l=0;l<dim;l++){
                //Find u_{k,l} for the current quadrature point
                du_dx[k][l] += D[local_dof_indices[elemDOF]]*fe_face_values.shape_grad(elemDOF,q)[l];
              }
            }
          }
          // calculate the shear stresses on face x1 = 1, storing in sigma_12 and sigma_13          
          for(unsigned int k = 0; k < dim; k++) {
            for(unsigned int l = 0; l < dim; l++) {
              sigma_13 += C(0, 2, k, l) * (du_dx[k][l] + du_dx[l][k]) / 2.;
              sigma_12 += C(0, 1, k, l) * (du_dx[k][l] + du_dx[l][k]) / 2.;
            }
          }
          torque += fe_face_values.JxW(q) * (sigma_13 * yarm - sigma_12 * zarm);          
        }
      }
    }
  }
  return torque;

}


/**************************************************************************************/
/*  This function directly solves "K_glocal * D = F_global" by implementing a matrix  */
/*  inverse - which is cost expensive, but of relatively high accuracy.               */
/**************************************************************************************/

template <int dim>
void FEM<dim>::solve(){

  SparseDirectUMFPACK  A;
  A.initialize(K);
  A.vmult (D, F); //D=K^{-1}*F

}


/**************************************************************************************/
/*  This function writes nodal solutions to "solution_Problem_Order.vtk" for          */
/*  assignment submission and results checking.                                       */
/**************************************************************************************/

template <int dim>
void FEM<dim>::output_results (){

  std::string filename;
  filename = "../sols/vtk/solution_p" + std::to_string(problem) 
        +"_order" + std::to_string(order)
        + ".vtk";

  //Write results to VTK file
  std::ofstream output1 (filename);
  DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector (D,
                            nodal_solution_names,
                            DataOut<dim>::type_dof_data,
                            nodal_data_component_interpretation);
  data_out.build_patches();
  data_out.write_vtk(output1);
  output1.close();
}



template class FEM<1>;
template class FEM<2>;
template class FEM<3>;