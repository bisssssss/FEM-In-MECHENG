//Include files
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include "FEM1.h"
#include "writeSolutions.h"
#include "test.h"

using namespace dealii;

//The main program, using the FEM class
int main (){
  try{

    //Define the number of elements as an input to "generate_mesh"
    std::vector<unsigned int> mesh = {3, 10, 100, 1000};
    std::vector<std::vector<double>> l2err(mesh.size());
    std::vector<std::vector<double>> h1err(mesh.size());

    test_all();

    for(unsigned int imesh = 0; imesh < mesh.size(); imesh++){
        
        unsigned int nmesh = mesh[imesh];
        
        deallog.depth_console (0);

        //Specify the basis function order: 1, 2, or 3
        for(unsigned int problem = 1; problem <= 4; problem++){

            //Specify the subproblem: 1, 2, 3, or 4
            for(unsigned int order = 1; order <= 3; order++){

                //Case description
                std::cout << "p" << problem << " " <<  nmesh;
                switch (order) {
                    case 1: std::cout << " linear: " << std::endl; break;
                    case 2: std::cout << " quadratic: " << std::endl; break;
                    case 3: std::cout << " cubic: " << std::endl; break;
                    default: std::cout << " ?undefined order?: " << std::endl; break;
                }

                // Initialize problem
                FEM<1> problemObject(order,problem);
                
                // Solving problem
                problemObject.generate_mesh(nmesh);
                problemObject.setup_system();
                problemObject.assemble_system();
                problemObject.solve();
                
                // Show L2 and H1 error
                std::cout << "   L2 error:                     "
                    << problemObject.l2norm_of_error() << std::endl;
                std::cout << "   H1 error:                     "
                    << problemObject.h1norm_of_error() << std::endl; std::cout << std::endl;
                l2err[imesh].push_back(problemObject.l2norm_of_error());
                h1err[imesh].push_back(problemObject.h1norm_of_error());

                //write output file in vtk format for visualization
                problemObject.output_results();
                
                //write solutions to h5 file
                char tag[21];
                sprintf(tag, "CA1_Order%d_Problem%d",order,problem);
                writeSolutionsToFileCA1(problemObject.D, problemObject.l2norm_of_error(), tag);

            } //for(order)

        } // for(problem)

    } // for(imesh)

    std::ofstream fname;
    fname.open ("../sols/plots/data/l2err.txt");
    for(const auto row : l2err) {
        for(const auto col : row)
            fname << col << " ";
        fname << "\n";
    }
    fname.close();

    std::ofstream fname2;
    fname2.open ("../sols/plots/data/h1err.txt");
    for(const auto row : h1err) {
        for(const auto col : row)
            fname2 << col << " ";
        fname2 << "\n";
    }
    fname2.close();
  } // try()

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