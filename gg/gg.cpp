// gg.cpp
// Grain Growth model for PPC project

#ifndef SOLIDIFICATION_UPDATE
#define SOLIDIFICATION_UPDATE
#include<iomanip>
#include<cmath>
#include<ctime>
#include"MMSP.hpp"
#include"gg.hpp"
#include "tessellate.hpp"

namespace MMSP
{

/* Grid contains three fields:
   0. Concentration (c)
   1. Phase (phi)
   2. State (s)
*/

void generate(int dim, const char* filename)
{
	const int L = 32; 
	const double deltaX = 0.025;
	const double init_c = 0.5;
	int number_of_fields = static_cast<int>(float(L * L) / (M_PI * 1 * 1)); /* average grain is a disk of radius XXX, XXX cannot be smaller than 0.1, or BGQ will abort.*/

#ifdef MPI_VERSION
		number_of_fields /= np;
#endif

	MMSP::grid<dim, unsigned long > grid = MMSP::grid<dim, unsigned long>(0, 0, L, 0, L);
	unsigned long timer = tessellate<dim, unsigned long>(grid, number_of_fields, 1);

	if (dim == 2) {
		GRID2D initGrid(3,0,L,0,L); //0 is concentration, 1 is phase, 2 is state.
		for (int d = 0; d < dim; d++)
			dx(initGrid,d) = deltaX;

		for (int i = 0; i < nodes(initGrid); i++) {
			initGrid(i)[2] = grid(i);
			if(phases[grid(i)] == 1) {
				initGrid(i)[1] = 1;
				initGrid(i)[0] = 0.7;
			}
			else {
				initGrid(i)[1] = 0;
				initGrid(i)[0] = 0.3;
			}
		}
		output(initGrid,filename);
	} else {
		std::cerr << "Grain Growth code is only implemented for 2D." << std::endl;
		MMSP::Abort(-1);
	}
}

template <typename T> T f(const T c, const T phi) { 
	// 0 = alpha, 1 = beta, for phi
	const double lmbda = 0.3;   
	const double   C_1 = 0.25;  
	const double   C_2 = 0.75;    
	const double   C_3 = 0.05; 
	const double   C_4 = 0.95;  
	const double     a = 0.5;
	return lmbda*(4*c-2*C_1-2*C_2) + a*(2*c-2*C_3)*(1-phi) + a*(2*c-2*C_4)*phi;
}

template <int dim, typename T> void update_conc(grid<dim,vector<T> >& oldGrid, int steps)
{
	int id = 0;
	int np = 1;
	static int iterations = 1;
	#ifdef MPI_VERSION
	id = MPI::COMM_WORLD.Get_rank();
	np = MPI::COMM_WORLD.Get_size();
	#endif

	ghostswap(oldGrid);

	grid<dim,vector<T> > newGrid(oldGrid);
	
	const double    dt = 5e-5;  
	const double   M_C = 0.1;
	const double   K_C = 0.001;
	
	std::cout.precision(2);

	ghostswap(oldGrid);

	int minus = 0;
	int plus = 0;
	for (int step = 0; step < steps; step++) {
		if (id == 0)
			print_progress(step, steps);

		grid<dim,vector<T> > c_storage(oldGrid); //0: laplacian_C, 1: laplacian_(f(C))

		for (int i = 0; i < nodes(oldGrid); i++) {
			vector<int> x = position(oldGrid,i);
			const T& center = oldGrid(x)[0];
			const T& center_E = f(center, oldGrid(x)[1]);
			x[0]++;
			const T& right = oldGrid(x)[0];
			const T& right_E = f(right, oldGrid(x)[1]);
			x[0] -= 2;
			const T& left = oldGrid(x)[0];
			const T& left_E = f(left, oldGrid(x)[1]);
			x[0]++;
			x[1]++;
			const T& up = oldGrid(x)[0];
			const T& up_E = f(up, oldGrid(x)[1]);
			x[1] -= 2;
			const T& down = oldGrid(x)[0];
			const T& down_E = f(down, oldGrid(x)[1]);
			x[1]++;
			c_storage(i)[0] = (left+right+up+down-4*center)/(dx(oldGrid, 0)*dx(oldGrid, 0));
			c_storage(i)[1] = (left_E+right_E+up_E+down_E-4*center_E)/(dx(oldGrid, 0)*dx(oldGrid, 0));	

		}

		// Sync parallel grids
		ghostswap(c_storage);
		ghostswap(oldGrid);

		for (int i = 0; i < nodes(oldGrid); i++) {
			vector<int> x = position(oldGrid,i);
			const vector<T>& oldX = oldGrid(x);
			const vector<T>&  csX = c_storage(x);
			      vector<T>& newX = newGrid(x);

			// Update concentration
			x[0]++;
			const T& right = c_storage(x)[0];
			x[0] -= 2;
			const T& left = c_storage(x)[0];
			x[0]++;
			x[1]++;
			const T& up = c_storage(x)[0];
			x[1] -= 2;
			const T& down = c_storage(x)[0];
			x[1]++;
			
			T Kg4C = K_C*(left+right+up+down-4*csX[0])/(dx(oldGrid, 0)*dx(oldGrid, 0));
			T g2dEdC = csX[1];
			T dCdt = M_C*(g2dEdC - Kg4C);
			newX[0] = oldX[0] + dt*dCdt;
			newX[1] = oldX[1];
			newX[2] = oldX[2];
		}

		swap(oldGrid,newGrid);
		ghostswap(oldGrid);
	}

	iterations++;
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
