//Include files
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "FEM4.h"
#include "writeSolutions.h"

using namespace dealii;

//The main program, using the FEM class
int main (){
  try{
    deallog.depth_console (0);
    std::vector<double> delta_tset;
    double delta_t = 1.;
    double alpha = 0.;

    //double delta_t;
    //for(double alpha = 0.; alpha < 1.1; alpha += 0.5) {
      //for(double dt = 0.1; dt < 15; dt*=1.5) {
        //delta_t = dt;

        const int dimension = 3;
        //double alpha = 1.; //Specify the Euler method (0 <= alpha <= 1)
        FEM<dimension> problemObject(alpha, delta_t);

        std::vector<unsigned int> num_of_elems(dimension);
        num_of_elems[0] = 20;
        num_of_elems[1] = 20;
        num_of_elems[2] = 1; //For example, a 10x10x1 mesh

        problemObject.generate_mesh(num_of_elems);
        problemObject.setup_system();
        problemObject.assemble_system();

        //Solve the steady state problem
        problemObject.solve_steady();

        //Apply solve the transient problems
        problemObject.solve_trans();

        //Print l2norms to screen
        for(unsigned int i=0; i<problemObject.l2norm_results.size(); i++){
          std::cout << problemObject.l2norm_results[i] << std::endl;
        }

        //write l2 norm vector to a txt file for post plotting
        /*
        std::ofstream outf;
        outf.open ("../sols/l2_0.txt");
        for(unsigned int i=0; i<problemObject.l2norm_results.size(); i++){
          outf << problemObject.l2norm_results[i] << "\n";
        }
        outf.close();
        */

        //write l2 norm vector to a txt file for stability analysis
        /*
        std::ofstream outf;
        outf.open ("../sols/stab/l2_a" + std::to_string(alpha) + "_t" + std::to_string(delta_t)+ ".txt");
        for(unsigned int i=0; i<problemObject.l2norm_results.size(); i++){
          outf << problemObject.l2norm_results[i] << "\n";
        }
        outf.close();
        */

        //write solutions to h5 file
        char tag[21];
        sprintf(tag, "CA4_Alpha%2.1f",alpha);
        writeSolutionsToFileCA4(problemObject.D_trans, problemObject.l2norm_results[30], tag);
      //}
    
    //}

  }
  catch (std::exception &exc){
    std::cerr << std::endl << std::endl
        << "----------------------------------------------------"
        << std::endl;
    std::cerr << "Exception on processing: " << std::endl
        << exc.what() << std::endl
        << "Aborting!" << std::endl
        << "----------------------------------------------------"
        << std::endl;

    return 1;
  }
  catch (...){
    std::cerr << std::endl << std::endl
        << "----------------------------------------------------"
        << std::endl;
    std::cerr << "Unknown exception!" << std::endl
        << "Aborting!" << std::endl
        << "----------------------------------------------------"
        << std::endl;
    return 1;
  }

  return 0;
}
