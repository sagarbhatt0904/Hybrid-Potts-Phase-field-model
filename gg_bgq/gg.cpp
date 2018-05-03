// gg.cpp
// Grain Growth model for PPC project

#include<iomanip>
#include<cmath>
#include<ctime>
#include"MMSP.hpp"
#include"gg.hpp"
#include <vector>
#include "tessellate.hpp"

#define L 1001

namespace MMSP
{


//function to compute dE/dC
template <typename T> T f(const T c, const T phi) {
	// 0 = alpha, 1 = beta, for phi
	const double lmbda = 0.3;
	const double   C_1 = 0.25;
	const double   C_2 = 0.75;
	const double   C_3 = 0.05;
	const double   C_4 = 0.95;
	const double     a = 0.5;
	return lmbda * (4 * c - 2 * C_1 - 2 * C_2) + a * (2 * c - 2 * C_3) * (1 - phi) + a * (2 * c - 2 * C_4) * phi;
}

//function to compute E_total, including effects from C, phi, and grain ID
template <typename T> std::vector<double> dE(const T c, const std::vector<double>& phases, const std::vector<double>& grains) {
	// 0 = alpha, 1 = beta, for phi
	const double lmbda = 0.3;
	const double   C_1 = 0.25;
	const double   C_2 = 0.75;
	const double   C_3 = 0.05;
	const double   C_4 = 0.95;
	const double     a = 0.5;
	std::vector<double> deltas;
	deltas.push_back(0);
	T original = lmbda * ((c - C_1) * (c - C_1) + (c - C_2) * (c - C_2)) + a * ((c - C_3) * (c - C_3)) * (1 - phases[0]) + a * ((c - C_4) * (c - C_4)) * phases[0];
	for (int k = 1; k < 5; k++) {
		original += (grains[k] != grains[0]);
	}
	for (int k = 1; k < 5; k++) {
		T changed = lmbda * ((c - C_1) * (c - C_1) + (c - C_2) * (c - C_2)) + a * ((c - C_3) * (c - C_3)) * (1 - phases[k]) + a * ((c - C_4) * (c - C_4)) * phases[k];
		for (int j = 1; j < 5; j++) {
			changed += (grains[k] != grains[j]);
		}
		deltas.push_back(changed - original);
	}
	return deltas;
}

//creates the grid, then simulates for [steps] steps
void update(int steps)
{
	//std::cout << "beginning of update\n";
	const double    dt = 5e-5;
	const double   M_C = 0.025;
	const double   K_C = 0.001;
	const double deltaX = 0.025;

	MMSP::grid<2, unsigned long > grid = MMSP::grid<2, unsigned long>(1, 0, L, 0, L);

	int number_of_fields = static_cast<int>(float(L * L) / (M_PI * 1 * 1)); /* average grain is a disk of radius XXX, XXX cannot be smaller than 0.1, or BGQ will abort.*/

#ifdef MPI_VERSION
	int np = MPI::COMM_WORLD.Get_size();
	number_of_fields /= np;
#endif
	//std::cout << "before tessellate\n";
	//unsigned long timer = tessellate<2, unsigned long>(grid, number_of_fields, 1);
	//std::cout << "after tessellate\n";
	GRID2D oldGrid(3, 0, L, 0, L); //0 is concentration, 1 is phase, 2 is state.
	for (int d = 0; d < 2; d++)
		dx(oldGrid, d) = deltaX;

	for (int i = 0; i < nodes(oldGrid); i++) {
		int phase = rand() % 2;
		int id = rand();
		oldGrid(i)[2] = id;
		if (phase == 1) {
			oldGrid(i)[1] = 1;
			oldGrid(i)[0] = 0.7;
		}
		else {
			oldGrid(i)[1] = 0;
			oldGrid(i)[0] = 0.3;
		}
	}
	
/* Grid contains three fields:
   0. Concentration (c)
   1. Phase (phi)
   2. State (s)
*/


	int id = 0;

	static int iterations = 1;
#ifdef MPI_VERSION
	id = MPI::COMM_WORLD.Get_rank();

#endif

	ghostswap(oldGrid);

	GRID2D newGrid(3, 0, L, 0, L);

	std::cout.precision(2);

	ghostswap(oldGrid);

	int minus = 0;
	int plus = 0;
	for (int step = 0; step < steps; step++) {
		if (id == 0)
			print_progress(step, steps);

		GRID2D c_storage(2, 0, L, 0, L); //0: laplacian_C, 1: laplacian_(f(C))

		for (int i = 0; i < nodes(oldGrid); i++) {
			vector<int> x = position(oldGrid, i);
			const double& center = oldGrid(x)[0];
			const double& center_E = f(center, oldGrid(x)[1]);
			x[0]++;
			const double& right = oldGrid(x)[0];
			const double& right_E = f(right, oldGrid(x)[1]);
			x[0] -= 2;
			const double& left = oldGrid(x)[0];
			const double& left_E = f(left, oldGrid(x)[1]);
			x[0]++;
			x[1]++;
			const double& up = oldGrid(x)[0];
			const double& up_E = f(up, oldGrid(x)[1]);
			x[1] -= 2;
			const double& down = oldGrid(x)[0];
			const double& down_E = f(down, oldGrid(x)[1]);
			x[1]++;
			c_storage(i)[0] = (left + right + up + down - 4 * center) / (dx(oldGrid, 0) * dx(oldGrid, 0));
			c_storage(i)[1] = (left_E + right_E + up_E + down_E - 4 * center_E) / (dx(oldGrid, 0) * dx(oldGrid, 0));

		}

		// Sync parallel grids
		ghostswap(c_storage);
		ghostswap(oldGrid);

		for (int i = 0; i < nodes(oldGrid); i++) {
			vector<int> x = position(oldGrid, i);
			const vector<double>& oldX = oldGrid(x);
			const vector<double>&  csX = c_storage(x);
			vector<double>& newX = newGrid(x);

			// Update concentration
			x[0]++;
			const double& right = c_storage(x)[0];
			x[0] -= 2;
			const double& left = c_storage(x)[0];
			x[0]++;
			x[1]++;
			const double& up = c_storage(x)[0];
			x[1] -= 2;
			const double& down = c_storage(x)[0];
			x[1]++;

			double Kg4C = K_C * (left + right + up + down - 4 * csX[0]) / (dx(oldGrid, 0) * dx(oldGrid, 0));
			double g2dEdC = csX[1];
			double dCdt = M_C * (g2dEdC - Kg4C);
			newX[0] = oldX[0] + dt * dCdt;
		}

		for (int i = 0; i < nodes(oldGrid); i++) {
			newGrid(i)[1] = oldGrid(i)[1];
			newGrid(i)[2] = oldGrid(i)[2];
			if (i % 2 == step % 2) {
				std::vector<double> grains; //0 = center, 1-4 are neighbors
				std::vector<double> phases; //same
				vector<int> x = position(oldGrid, i);
				grains.push_back(oldGrid(x)[2]);
				phases.push_back(oldGrid(x)[1]);
				x[0]++;
				grains.push_back(oldGrid(x)[2]);
				phases.push_back(oldGrid(x)[1]);
				x[0] -= 2;
				grains.push_back(oldGrid(x)[2]);
				phases.push_back(oldGrid(x)[1]);
				x[0]++;
				x[1]++;
				grains.push_back(oldGrid(x)[2]);
				phases.push_back(oldGrid(x)[1]);
				x[1] -= 2;
				grains.push_back(oldGrid(x)[2]);
				phases.push_back(oldGrid(x)[1]);
				x[1]++;
				std::vector<double> deltas = dE(oldGrid(x)[0], phases, grains);
				int j = rand() % 4;
				if (deltas[j + 1] < 0) {
					newGrid(x)[1] = phases[j + 1];
					newGrid(x)[2] = grains[j + 1];
				}
				else if ((((double)rand()) / RAND_MAX) < exp(-deltas[j + 1])) {
					newGrid(x)[1] = phases[j + 1];
					newGrid(x)[2] = grains[j + 1];
				}
			}

		}

		swap(oldGrid, newGrid);
		ghostswap(oldGrid);
	}

	iterations++;
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
