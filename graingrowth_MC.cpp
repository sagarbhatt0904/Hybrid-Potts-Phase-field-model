#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <pthread.h>
#include "rdtsc.h"
#include"graingrowth_MC.hpp"
#include"MMSP.hpp"
#include"tessellate.hpp"
#include"output.cpp"
#include "gg.cpp"

// grid point dimension
int L = 64;

#ifndef SILENT
void print_progress(const int step, const int steps, const int iterations) {
	char* timestring;
	static unsigned long tstart;
	struct tm* timeinfo;

	if (step == 0) {
		tstart = time(NULL);
		std::time_t rawtime;
		std::time( &rawtime );
		timeinfo = std::localtime( &rawtime );
		timestring = std::asctime(timeinfo);
		timestring[std::strlen(timestring) - 1] = '\0';
		std::cout << "Pass " << std::setw(3) << std::right << iterations << ": " << timestring << " [" << std::flush;
	} else if (step == steps) {
		unsigned long deltat = time(NULL) - tstart;
		std::cout << "•] "
		          << std::setw(2) << std::right << deltat / 3600 << "h:"
		          << std::setw(2) << std::right << (deltat % 3600) / 60 << "m:"
		          << std::setw(2) << std::right << deltat % 60 << "s"
		          << " (File " << std::setw(5) << std::right << iterations*steps << ")." << std::endl;
	} else if ((20 * step) % steps == 0) std::cout << "• " << std::flush;
}
#endif


namespace MMSP
{
template <int dim> bool OutsideDomainCheck(MMSP::grid<dim, unsigned long>& grid, vector<int>* x);

template <int dim>
unsigned long generate(MMSP::grid<dim, unsigned long >*& grid, int seeds, int nthreads, const char* file2)
{
#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr << "Error: MPI is required for CCNI." << std::endl;
	exit(1);
#endif
#ifdef MPI_VERSION
	int np = MPI::COMM_WORLD.Get_size();
#endif

	unsigned long timer = 0;
	if (dim == 2) {
		int number_of_fields(seeds);
		if (number_of_fields == 0) number_of_fields = static_cast<int>(float(L * L) / (M_PI * 1 * 1)); /* average grain is a disk of radius XXX, XXX cannot be smaller than 0.1, or BGQ will abort.*/
#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
#endif
		grid = new MMSP::grid<dim, unsigned long>(3, 0, L, 0, L);

#ifdef MPI_VERSION
		number_of_fields /= np;
#endif

#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr << "Error: CCNI requires MPI." << std::endl;
		std::exit(1);
#endif
		timer = tessellate<dim, unsigned long>(*grid, number_of_fields, nthreads);
		vector<int> r(dim, 0);
		//if (id == 0) seed_points << "x,y\n";
		generate_conc( dim, grid, file2);





#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
#endif
	}

	return timer;
}



int LargeNearestInteger(int a, int b) {
	if (a % b == 0) return a / b;
	else return a / b + 1;
}

template <int dim> struct flip_index {
	MMSP::grid<dim, unsigned long>* grid;
	MMSP::grid<dim, MMSP::vector<unsigned long>>&  oldGrid;
	int num_of_cells_in_thread;
	int sublattice;
	int num_of_points_to_flip;
	int cell_coord[dim];
	int lattice_cells_each_dimension[dim];
	double Pdenominator;
};

template <int dim> void* flip_index_helper( void* s )
{
	srand(time(NULL)); /* seed random number generator */
	flip_index<dim>* ss = static_cast<flip_index<dim>*>(s);
	//  double kT=0.0;
	vector<int> x (dim, 0);
	int first_cell_start_coordinates[dim];
	for (int kk = 0; kk < dim; kk++) first_cell_start_coordinates[kk] = x0(*(ss->grid), kk);
	for (int i = 0; i < dim; i++) {
		if (x0(*(ss->grid), i) % 2 != 0) first_cell_start_coordinates[i]--;
	}
	int cell_coords_selected[dim];
	for (int hh = 0; hh < ss->num_of_points_to_flip; hh++) {
		// choose a random cell to flip
		int cell_numbering_in_thread = rand() % (ss->num_of_cells_in_thread); //choose a cell to flip, from 0 to num_of_cells_in_thread-1
		if (dim == 2) {
			cell_coords_selected[dim - 1] = ((ss->cell_coord)[dim - 1] + cell_numbering_in_thread) % (ss->lattice_cells_each_dimension)[dim - 1]; //1-indexed
			cell_coords_selected[0] = (ss->cell_coord)[0] + (((ss->cell_coord)[dim - 1] + cell_numbering_in_thread) / (ss->lattice_cells_each_dimension)[dim - 1]);
		}

		for (int i = 0; i < dim; i++) {
			x[i] = first_cell_start_coordinates[i] + 2 * cell_coords_selected[i];
		}
		if (dim == 2) {
			switch (ss->sublattice) {
			case 0: break; // 0,0
			case 1: x[1]++; break; //0,1
			case 2: x[0]++; break; //1,0
			case 3: x[0]++; x[1]++; break; //1,1
			}
		}


		bool site_out_of_domain = false;
		for (int i = 0; i < dim; i++) {
			if (x[i] < x0(*(ss->grid), i) || x[i] > x1(*(ss->grid), i)) {
//      if(x[i]<x0(*(ss->grid), i) || x[i]>x1(*(ss->grid), i)-1){
				site_out_of_domain = true;
				break;//break from the for int i loop
			}
		}
		if (site_out_of_domain == true) {
			hh--;
			continue; //continue the int hh loop
		}

		int rank = MPI::COMM_WORLD.Get_rank();

		unsigned long spin1 = (*(ss->grid))(x);
		unsigned long phase1 = oldGrid(x)[1];
		// determine neighboring spins
		vector<int> r(dim, 0);
		std::vector<unsigned long> neighbors;
		neighbors.clear();
		std::vector<unsigned long> neighborsPhase;
		neighborsPhase.clear();
		int number_of_same_neighours = 0;
		if (dim == 2) {
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					if (!(i == 0 && j == 0)) {
						r[0] = x[0] + i;
						r[1] = x[1] + j;
						unsigned long spin = (*(ss->grid))(r);
						unsigned long phase = oldGrid(r)[1];
						neighbors.push_back(spin);
						neighborsPhase.push_back(phase);
						if (spin == spin1)
							number_of_same_neighours++;
					}
				}
			}
		}


		//check if inside a grain
		if (number_of_same_neighours == neighbors.size()) { //inside a grain
			continue;//continue for
		}
		//choose a random neighbor spin
		unsigned long spin2 = neighbors[rand() % neighbors.size()];

		if (spin1 != spin2) {
			// compute energy change
			double dE = 0.0;
			if (dim == 2) {
				for (int i = -1; i <= 1; i++) {
					for (int j = -1; j <= 1; j++) {
						if (!(i == 0 && j == 0)) {
							r[0] = x[0] + i;
							r[1] = x[1] + j;
							unsigned long spin = (*(ss->grid))(r);
							dE += 1.0 / 2 * ((spin != spin2) - (spin != spin1));
						}// if(!(i==0 && j==0))
					}
				}
			}

			if (dE <= 0.0) {
				(*(ss->grid))(x) = spin2;
			}

		}
	}
	pthread_exit(0);
	return NULL;
}


template <int dim> bool OutsideDomainCheck(MMSP::grid<dim, unsigned long>& grid, vector<int>* x) {
	bool outside_domain = false;
	for (int i = 0; i < dim; i++) {
		if ((*x)[i] < x0(grid, i) || (*x)[i] > x1(grid, i) - 1) {
			outside_domain = true;
			break;
		}
	}
	return outside_domain;
}

template <int dim, typename T> unsigned long update(MMSP::grid<dim, unsigned long>& grid,  MMSP::grid<dim, vector<T> >& oldGrid, int steps, int nthreads)
{

#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	std::cerr << "Error: MPI is required for CCNI." << std::endl;
	exit(1);
#endif
	int rank = 0;
	unsigned int np = 0;
//	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	np = MPI::COMM_WORLD.Get_size();
	MPI::COMM_WORLD.Barrier();
//	#endif

	unsigned long update_timer = 0;

	pthread_t* p_threads = new pthread_t[nthreads];
	flip_index<dim>* mat_para = new flip_index<dim> [nthreads];
	pthread_attr_t attr;
	pthread_attr_init (&attr);

#ifndef SILENT
	static int iterations = 1;
	if (rank == 0) print_progress(0, steps, iterations);
#endif

	update_conc(oldGrid, steps);


	/*-----------------------------------------------*/
	/*---------------generate cells------------------*/
	/*-----------------------------------------------*/
	int dimension_length = 0, number_of_lattice_cells = 1;
	int lattice_cells_each_dimension[dim];
	for (int i = 0; i < dim; i++) {
		dimension_length = x1(grid, i) - x0(grid, i);
//    dimension_length = x1(grid, i)-1-x0(grid, i);
		if (x0(grid, 0) % 2 == 0)
			lattice_cells_each_dimension[i] = dimension_length / 2 + 1;
		else
			lattice_cells_each_dimension[i] = 1 + LargeNearestInteger(dimension_length, 2);
		number_of_lattice_cells *= lattice_cells_each_dimension[i];
	}
#ifndef SILENT
	if (rank == 0) {
		std::cout << "lattice_cells_each_dimension are:  ";
		for (int i = 0; i < dim; i++)
			std::cout << lattice_cells_each_dimension[i] << "  ";
		std::cout << "\n";
	}
#endif
//----------assign cells for each pthreads
	int num_of_cells_in_thread = number_of_lattice_cells / nthreads;
	//check if num of the pthread is too large, if so, reduce it.
	if (num_of_cells_in_thread < 1) {
		std::cerr << "ERROR: number of pthread is too large, please reduce it to a value <= " << number_of_lattice_cells << std::endl;
		exit(0);
	}

	vector<int> x (dim, 0);
	vector<int> x_prim (dim, 0);
	int coordinates_of_cell[dim];
	int initial_coordinates[dim];
	int **cell_coord = new int*[nthreads];//record the start coordinates of each pthread domain.
	for (int i = 0; i < nthreads; i++) {
		cell_coord[i] = new int[dim];
		for (int j = 0; j < dim; j++) {
			cell_coord[i][j] = 0;
		}
	}

	int **num_of_grids_to_flip = new int*[nthreads];
	for (int i = 0; i < nthreads; i++) {
		num_of_grids_to_flip[i] = new int[( static_cast<int>(pow(2, dim)) )];
		for (int j = 0; j < pow(2, dim); j++) {
			num_of_grids_to_flip[i][j] = 0;
		}
	}

	for (int k = 0; k < dim; k++)
		initial_coordinates[k] = x0(grid, k);
	for (int i = 0; i < dim; i++) {
		if (x0(grid, i) % 2 != 0)
			initial_coordinates[i]--;
	}

	for (int i = 0; i < nthreads; i++) {
		int cell_numbering = num_of_cells_in_thread * i; //0-indexed, celling_numbering is the start cell numbering
		if (dim == 2) {
			cell_coord[i][dim - 1] = cell_numbering % lattice_cells_each_dimension[dim - 1]; //0-indexed
			cell_coord[i][0] = (cell_numbering / lattice_cells_each_dimension[dim - 1]);
			if (cell_coord[i][0] >= lattice_cells_each_dimension[0]) {
				std::cerr << "the cell coordinates is wrong!" << std::endl;
				exit(1);
			}
		}


		mat_para[i].grid = &grid;
		mat_para[i].oldGrid = oldGrid;
		if (i == (nthreads - 1))
			mat_para[i].num_of_cells_in_thread = number_of_lattice_cells - num_of_cells_in_thread * (nthreads - 1);
		else
			mat_para[i].num_of_cells_in_thread = num_of_cells_in_thread;

#ifndef SILENT
		if (rank == 0) std::cout << "num_of_cells_in_thread is " << mat_para[i].num_of_cells_in_thread << " in thread " << i << "\n";
#endif

		for (int k = 0; k < dim; k++)
			mat_para[i].lattice_cells_each_dimension[k] = lattice_cells_each_dimension[k];

		for (int j = 0; j < mat_para[i].num_of_cells_in_thread; j++) {
			int start_cell_numbering = num_of_cells_in_thread * i;
			if (dim == 2) {
				coordinates_of_cell[dim - 1] = (start_cell_numbering + j) % lattice_cells_each_dimension[dim - 1]; //0-indexed
				coordinates_of_cell[0] = (start_cell_numbering + j) / lattice_cells_each_dimension[dim - 1];
			}

			for (int ii = 0; ii < dim; ii++) {
				x[ii] = initial_coordinates[ii] + 2 * coordinates_of_cell[ii];
			}

			if (dim == 2) {
				x_prim = x;
				if (!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][0] += 1;

				x_prim = x;
				x_prim[1] = x[1] + 1; //0,1
				if (!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][1] += 1;

				x_prim = x;
				x_prim[0] = x[0] + 1; //1,0
				if (!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][2] += 1;

				x_prim = x;
				x_prim[0] = x[0] + 1;
				x_prim[1] = x[1] + 1; //1,1
				if (!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][3] += 1;
			}

		}// for int j
	}//for int i

	for (int step = 0; step < steps; step++) {
		unsigned long start = rdtsc();
		int num_of_sublattices = 0;
		if (dim == 2) num_of_sublattices = 4;
		for (int sublattice = 0; sublattice < num_of_sublattices; sublattice++) {
			for (int i = 0; i != nthreads ; i++) {
				mat_para[i].sublattice = sublattice;
				mat_para[i].num_of_points_to_flip = num_of_grids_to_flip[i][sublattice];
				for (int k = 0; k < dim; k++) mat_para[i].cell_coord[k] = cell_coord[i][k];
				pthread_create(&p_threads[i], &attr, flip_index_helper<dim>, (void*) &mat_para[i] );
			}//loop over threads

			for (int ii = 0; ii != nthreads ; ii++)
				pthread_join(p_threads[ii], NULL);

#ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
#endif
			if (rank == 0)
			{
				for (int i = 0; i < nodes(oldGrid); i++) {
					vector<int> x = position(oldGrid, i);
					oldGrid(x)[2] = (*(grid))(i);
				}

//			ghostswap(grid, sublattice); // once looped over a "color", ghostswap.
				ghostswap(grid); // once looped over a "color", ghostswap.
#ifdef MPI_VERSION
				MPI::COMM_WORLD.Barrier();
#endif

			}//loop over color
#ifndef SILENT
			if (rank == 0) print_progress(step + 1, steps, iterations);
#endif
			update_timer += rdtsc() - start;
		}//loop over step
#ifndef SILENT
		++iterations;
#endif

		for (int i = 0; i < nthreads; i++) {
			delete [] num_of_grids_to_flip[i];
			num_of_grids_to_flip[i] = NULL;
			delete [] cell_coord[i];
			cell_coord[i] = NULL;
		}
		delete [] num_of_grids_to_flip;
		num_of_grids_to_flip = NULL;
		delete [] cell_coord;
		cell_coord = NULL;

		delete [] p_threads;
		p_threads = NULL;
		delete [] mat_para;
		mat_para = NULL;
		unsigned long total_update_time = update_timer;
#ifdef MPI_VERSION
		MPI::COMM_WORLD.Allreduce(&update_timer, &total_update_time, 1, MPI_UNSIGNED_LONG, MPI_SUM);
#endif
		return total_update_time / np; // average update time
	}

}


#endif

#include"MMSP.main.hpp"

