//Include files
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include "FEM2a.h"
#include "writeSolutions.h"
#include "test2a.h"

using namespace dealii;

//The main program, using the FEM class
int main (){
  try{
    deallog.depth_console (0);

    const int dimension = 2;

    test_all();

    //Define basis order: 1 or 2 or 3 or higher
    int order = 1;

		//Define problem number: 1 or 2
		unsigned int problem = 1;

    //Define the number of elements in the mesh
    //Here we use a 15 x 40 element mesh in 2D
    std::vector<unsigned int> num_of_elems(dimension);
    num_of_elems[0] = 15;
    num_of_elems[1] = 40;

    //Initialize and run simulation
    FEM<dimension> problemObject(order, problem);
    problemObject.generate_mesh(num_of_elems);
    problemObject.setup_system();
    problemObject.assemble_system();
    problemObject.solve();
    problemObject.output_results();

		if(problem == 2){
			std::cout << problemObject.l2norm_of_error() << std::endl;
		}

    //write solutions to h5 file
    char tag[21];
    sprintf(tag, "CA2a_Problem%d",problem);
		if(problem == 1){
	    writeSolutionsToFileCA2_1(problemObject.D, tag);
		}
		else if(problem == 2){
	    writeSolutionsToFileCA2_2(problemObject.D, problemObject.l2norm_of_error(), tag);
		}
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
