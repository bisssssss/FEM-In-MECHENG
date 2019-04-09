//Include files
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include "FEM2b.h"
#include "writeSolutions.h"
#include "test2b.h"

using namespace dealii;

//The main program, using the FEM class
int main (){
  try{
    deallog.depth_console (0);
		
	const int dimension = 3;

	//test_all();

	//Define basis order: 1 or 2 or 3 or higher
    int order = 2;

	//Define the number of elements in the mesh
    //Here we use a 8 x 16 x 4 element mesh in 2D
	std::vector<unsigned int> num_of_elems(dimension);
	num_of_elems[0] = 10;
	num_of_elems[1] = 10;
	num_of_elems[2] = 10;

	//Initialize and run simulation
	FEM<dimension> problemObject(order);
	problemObject.generate_mesh(num_of_elems);
	problemObject.setup_system();
	problemObject.assemble_system();
	problemObject.solve();
	problemObject.output_results();
    
    //write solutions to h5 file
    char tag[21];
    sprintf(tag, "CA2b");
    writeSolutionsToFileCA2(problemObject.D, tag);

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
