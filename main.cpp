/***************************************************************************/
/* Code for Parallel Computing Project**************************************/
/* Scott Peters, Sagar Bhatt   *********************************************/
/*CSCI 6360: Parallel Computing*********************************************/
/***************************************************************************/


/***************************************************************************/
// After Copmiling, run <executable> --help for instructions
// Or,
// Run  <executable> <input_file> #number_of_iterations

/***************************************************************************/

#ifndef MMSP_MAIN
#define MMSP_MAIN
#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<cctype>
#include<time.h>
#include <mpi.h>
#include "gg.cpp"
#include"rdtsc.h"

// #define BGQ 1 // when running BG/Q, comment out when running on kratos
#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
double processor_frequency = 1600000000.0;
unsigned long long start_cycles = 0;
unsigned long long end_cycles = 0;
#else
#define GetTimeBase MPI_Wtime
double start_cycles = 0.0, end_cycles = 0.0;
#endif
double time_in_secs1 = 0.0;

int main(int argc, char* argv[]) {
	MMSP::Init(argc, argv);

	// check argument list
	if (argc < 2) {
		std::cout << PROGRAM << ": bad argument list.  Use\n\n";
		std::cout << "    " << PROGRAM << " --help\n\n";
		std::cout << "to generate help message.\n\n";
		MMSP::Abort(-1);
	}

	int rank = 0;
#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
#endif
	srand(time(NULL) + rank);
	/***************************************************************************/
	// print help message and exit
	/***************************************************************************/
	if (std::string(argv[1]) == std::string("--help")) {
		std::cout << PROGRAM << ": " << MESSAGE << "\n\n";
		std::cout << "Valid command lines have the form:\n\n";
		std::cout << "    " << PROGRAM << " ";
		std::cout << "[--help] [--example dimension [outfile]] [infile [outfile] steps [increment]]\n\n";
		std::cout << "A few examples of using the command line follow.\n\n";
		std::cout << "The command\n\n";
		std::cout << "    " << PROGRAM << " --help\n\n";
		std::cout << "generates this help message and exits.  ";
		std::cout << "The \"--example\" option can be used to gen-\nerate a relevant test grid, e.g.\n\n";
		std::cout << "    " << PROGRAM << " --example 3\n\n";
		std::cout << "generates an example test problem on a grid of dimension 3 and writes this to the \n";
		std::cout << "file named \"example\", while\n\n";
		std::cout << "    " << PROGRAM << " --example 2 start\n\n";
		std::cout << "generates an example test problem on a grid of dimension 2 and writes this to the \n";
		std::cout << "file named \"start\".\n\n";
		std::cout << "    " << PROGRAM << " start 1000\n\n";
		std::cout << "reads the grid contained within \"start\" and runs a simulation for 1000 time steps.\n";
		std::cout << "The final grid is written to a file named \"start.1000\".\n\n";
		std::cout << "    " << PROGRAM << " start final 1000\n\n";
		std::cout << "reads the grid contained within \"start\" and runs a simulation for 1000 time steps.\n";
		std::cout << "The final grid is written to a file named \"final.1000\".\n\n";
		std::cout << "    " << PROGRAM << " start 1000 100\n\n";
		std::cout << "reads the grid contained within \"start\" and runs a simulation for 1000 time steps.\n";
		std::cout << "The grid is then written to a file every 100 time steps.  ";
		std::cout << "The resulting files are \nnamed \"start.0100\", \"start.0200\", ... \"start.1000\".\n\n";
		std::cout << "    " << PROGRAM << " start final 1000 100\n\n";
		std::cout << "reads the grid contained within \"start\" and runs a simulation for 1000 time steps.\n";
		std::cout << "The grid is then written to a file every 100 time steps.  ";
		std::cout << "The resulting files are \nnamed \"final.0100\", \"final.0200\", ... \"final.1000\".\n\n";
		exit(0);
	}

	/***************************************************************************/
	// generate example grid
	/***************************************************************************/
	else if (std::string(argv[1]) == std::string("--example")) {
		// check argument list
		if (argc<3 or argc>4) {
			std::cout << PROGRAM << ": bad argument list.  Use\n\n";
			std::cout << "    " << PROGRAM << " --help\n\n";
			std::cout << "to generate help message.\n\n";
			MMSP::Abort(-1);
		}

		/***************************************************************************/
		// check problem dimension
		/***************************************************************************/
		if (std::string(argv[2]).find_first_not_of("0123456789") != std::string::npos) {
			std::cout << PROGRAM << ": example grid must have integral dimension.  Use\n\n";
			std::cout << "    " << PROGRAM << " --help\n\n";
			std::cout << "to generate help message.\n\n";
			MMSP::Abort(-1);
		}

		int dim = atoi(argv[2]);

		/***************************************************************************/
		// set output file name
		/***************************************************************************/
		std::string outfile;
		if (argc < 4) outfile = "example";
		else outfile = argv[3];

		char* filename = new char[outfile.length() + 1];
		for (unsigned int i = 0; i < outfile.length(); i++)
			filename[i] = outfile[i];
		filename[outfile.length()] = '\0';

		/***************************************************************************/
		// generate test problem
		/***************************************************************************/
		MMSP::generate(dim, filename);

		delete [] filename;
	}

	else {
		// bad argument list
		if (argc<3 or argc>5) {
			std::cout << PROGRAM << ": bad argument list.  Use\n\n";
			std::cout << "    " << PROGRAM << " --help\n\n";
			std::cout << "to generate help message.\n\n";
			MMSP::Abort(-1);
		}

		int steps;
		int increment;
		std::string outfile;

		if (std::string(argv[2]).find_first_not_of("0123456789") == std::string::npos) {
			// set output file name
			outfile = argv[1];

			// must have integral number of time steps
			if (std::string(argv[2]).find_first_not_of("0123456789") != std::string::npos) {
				std::cout << PROGRAM << ": number of time steps must have integral value.  Use\n\n";
				std::cout << "    " << PROGRAM << " --help\n\n";
				std::cout << "to generate help message.\n\n";
				MMSP::Abort(-1);
			}

			steps = atoi(argv[2]);
			increment = steps;

			if (argc > 3) {
				// must have integral output increment
				if (std::string(argv[3]).find_first_not_of("0123456789") != std::string::npos) {
					std::cout << PROGRAM << ": output increment must have integral value.  Use\n\n";
					std::cout << "    " << PROGRAM << " --help\n\n";
					std::cout << "to generate help message.\n\n";
					MMSP::Abort(-1);
				}

				increment = atoi(argv[3]);

				// output increment must be smaller than number of steps
				if (increment > steps) {
					std::cout << PROGRAM << ": output increment must be smaller than number of time steps.  Use\n\n";
					std::cout << "    " << PROGRAM << " --help\n\n";
					std::cout << "to generate help message.\n\n";
					MMSP::Abort(-1);
				}
			}
		}

		else {
			// set output file name
			outfile = argv[2];

			// set number of time steps
			if (std::string(argv[3]).find_first_not_of("0123456789") != std::string::npos) {
				// must have integral number of time steps
				std::cout << PROGRAM << ": number of time steps must have integral value.  Use\n\n";
				std::cout << "    " << PROGRAM << " --help\n\n";
				std::cout << "to generate help message.\n\n";
				MMSP::Abort(-1);
			}

			steps = atoi(argv[3]);
			increment = steps;

			if (argc > 4) {
				// must have integral output increment
				if (std::string(argv[4]).find_first_not_of("0123456789") != std::string::npos) {
					std::cout << PROGRAM << ": output increment must have integral value.  Use\n\n";
					std::cout << "    " << PROGRAM << " --help\n\n";
					std::cout << "to generate help message.\n\n";
					MMSP::Abort(-1);
				}

				increment = atoi(argv[4]);

				// output increment must be smaller than number of steps
				if (increment > steps) {
					std::cout << PROGRAM << ": output increment must be smaller than number of time steps.  Use\n\n";
					std::cout << "    " << PROGRAM << " --help\n\n";
					std::cout << "to generate help message.\n\n";
					MMSP::Abort(-1);
				}
			}
		}
	
		
		// file open error check

		std::ifstream input(argv[1]);
		if (!input) {
			std::cerr << "File input error: could not open " << argv[1] << ".\n\n";
			MMSP::Abort(-1);
		}

		// read data type
		std::string type;
		getline(input, type, '\n');

		// grid type error check
		if (type.substr(0, 4) != "grid") {
			std::cerr << "File input error: file does not contain grid data." << std::endl;
			MMSP::Abort(-1);
		}

		// read grid dimension
		int dim;
		input >> dim;

		input.close();
		/***************************************************************************/
		//Run Simulation
		/***************************************************************************/

		unsigned int rank = 0;
		unsigned int np = 0;
#ifdef MPI_VERSION
		rank = MPI::COMM_WORLD.Get_rank();
		np = MPI::COMM_WORLD.Get_size();
#endif
		if (dim == 2) {
			// construct grid object
			GRID2D grid(argv[1]);
			// Calculate Start time
			start_cycles = GetTimeBase();   //  Time start
			/***************************************************************************/
			// perform computation
			/***************************************************************************/
			MMSP::update(grid, steps, increment);


			end_cycles = GetTimeBase();   // Time end
#ifdef BGQ

			time_in_secs1 = ((double)(end_cycles - start_cycles)) / processor_frequency;
#else

			time_in_secs1 = ((double)(end_cycles - start_cycles));
#endif

			if (rank == 0)
			{
				/**********************************************************
				  Writing time to file for Plotting
				**********************************************************/
				FILE *ptrFileOut;
				ptrFileOut = fopen("time.dat", "a");
				fprintf( ptrFileOut, "%d %f \n", np, time_in_secs1);
				fclose(ptrFileOut);
			}
			MMSP::output(grid, "bgfinal");

		}
	}

	MMSP::Finalize();
}

#endif
