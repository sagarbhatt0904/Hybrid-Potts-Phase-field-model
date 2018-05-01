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


//----------Ti6Al4V Alloy

double lambda = 5.0e-3; 
double L_initial = 40.0e-3; 
double L0 = 40.0e-3; 
double K1 = 0.01; //Take a look at this  number as well
double n1 = 0.55;
double m = 1.0/n1;
double Q = 312000; 
double n = 2; 
double K_ = 1.11e+9; 
double R = 8.314; 




// grid point dimension
int dim_x =50; 
int dim_y = 50; 
int dim_z = 5; 
std::vector<std::vector<double> > u_t(dim_x,std::vector<double>(dim_y, 0));


#ifndef SILENT
void print_progress(const int step, const int steps, const int iterations) {
  char* timestring;
  static unsigned long tstart;
  struct tm* timeinfo;

  if (step==0) {
    tstart = time(NULL);
    std::time_t rawtime;
    std::time( &rawtime );
    timeinfo = std::localtime( &rawtime );
    timestring = std::asctime(timeinfo);
    timestring[std::strlen(timestring)-1] = '\0';
    std::cout<<"Pass "<<std::setw(3)<<std::right<<iterations<<": "<<timestring<<" ["<<std::flush;
  } else if (step==steps) {
    unsigned long deltat = time(NULL)-tstart;
    std::cout << "•] "
              <<std::setw(2)<<std::right<<deltat/3600<<"h:"
              <<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
              <<std::setw(2)<<std::right<<deltat%60<<"s"
              <<" (File "<<std::setw(5)<<std::right<<iterations*steps<<")."<<std::endl;
  } else if ((20 * step) % steps == 0) std::cout<<"• "<<std::flush;
}
#endif


namespace MMSP
{
template <int dim> bool OutsideDomainCheck(MMSP::grid<dim, unsigned long>& grid, vector<int>* x);

/*template <int dim> void UpdateLocalTmp(MMSP::grid<dim, unsigned long>& grid, long double physical_time, double* temp){
   vector<int> coords (dim,0);
   if(dim==2){
    
       for(int codx=x0(grid, 0); codx < x1(grid, 0); codx++) 
         for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
           coords[0] = codx;
           coords[1] = cody;
          grid.AccessToTmp(coords) = 1361.0;
           }

   }
   else if(dim==3){
     for(int codx=x0(grid, 0); codx < x1(grid, 0); codx++)  
         for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++) 
           for(int codz=x0(grid, 2); codz < x1(grid, 2); codz++){
             coords[0] = codx;
             coords[1] = cody;
             coords[2] = codz;
             grid.AccessToTmp(coords) = 473 + 273.0*sin(3.14/1000*codx + 3.14/(1.0e5)*physical_time);
           }
   }
}*/

template <int dim>
unsigned long generate(MMSP::grid<dim,unsigned long >*& grid, int seeds, int nthreads)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	#ifdef MPI_VERSION
	int np = MPI::COMM_WORLD.Get_size();
	#endif

	unsigned long timer=0;
	if (dim == 2) {
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(dim_x*dim_y)/(M_PI*1*1)); /* average grain is a disk of radius XXX
, XXX cannot be smaller than 0.1, or BGQ will abort.*/
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields; 
		#endif
		grid = new MMSP::grid<dim,unsigned long>(0, 0, dim_x, 0, dim_y);

		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
		timer = tessellate<dim,unsigned long>(*grid, number_of_fields, nthreads);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	} else if (dim == 3) {
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(dim_x*dim_y*dim_z)/(4./3*M_PI*0.5*0.5*0.5)); // Average grain is a sphere of radius 10 voxels
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
		#endif
		grid = new MMSP::grid<dim,unsigned long>(0, 0, dim_x, 0, dim_y, 0, dim_z);

		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

//		timer = tessellate<dim,unsigned long >(*grid, number_of_fields, nthreads);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	}
/*------------------Initial tmc----------------------*/
	double tmc_initial = 0;//pow((L_initial/lambda-1)/K1,1.0/n1);
  vector<int> coords (dim,0);
  if(dim==2){
    for(int codx=x0(*grid, 0); codx < x1(*grid, 0); codx++) 
      for(int cody=x0(*grid, 1); cody < x1(*grid, 1); cody++){
        coords[0] = codx;
        coords[1] = cody;
        // (*grid)(coords) = coords[0]*dim_y + coords[1] + 1;// grain id start from 1
        (*grid).AccessToTmc(coords) = tmc_initial;
      }
  }
  else if(dim==3){
    for(int codx=x0(*grid, 0); codx < x1(*grid, 0); codx++) 
      for(int cody=x0(*grid, 1); cody < x1(*grid, 1); cody++) 
        for(int codz=x0(*grid, 2); codz < x1(*grid, 2); codz++){
          coords[0] = codx;
          coords[1] = cody;
          coords[2] = codz;
          (*grid)(coords) = coords[0]*dim_y*dim_z + coords[1]*dim_z + coords[2] + 1;// grain id start from 1
          (*grid).AccessToTmc(coords) = tmc_initial;
        }
  }
  // double temp[2] = {673,3};
  //  UpdateLocalTmp((*grid), 0.0, &temp[0]);  
/*---------------------------------------------------*/
	return timer;
}

unsigned long generate(int dim, char* filename, int seeds, int nthreads)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	unsigned long timer = 0;
/*------------------Initial tmc----------------------*/
	double tmc_initial = 0;//pow((L_initial/lambda-1)/K1,1.0/n1);
  vector<int> coords (dim,0);
/*---------------------------------------------------*/

	if (dim == 2) {
		MMSP::grid<2,unsigned long>* grid2=NULL;
		timer = generate<2>(grid2,seeds,nthreads);
		assert(grid2!=NULL);
		#ifdef BGQ
		output_bgq(*grid2, filename);
		#else
		output(*grid2, filename);
		#endif
		#ifndef SILENT
		if (rank==0) std::cout<<"Wrote initial file to "<<filename<<"."<<std::endl;
		#endif
/*------------------Initial tmc----------------------*/
    for(int codx=x0(*grid2, 0); codx < x1(*grid2, 0); codx++) 
      for(int cody=x0(*grid2, 1); cody < x1(*grid2, 1); cody++){
        coords[0] = codx;
        coords[1] = cody;
        (*grid2).AccessToTmc(coords) = tmc_initial;
//        (*grid2).AccessToTmp(coords) = 1.0e6; //set the initial temp to be large enough such that initial physical time = 0.
      }
/*---------------------------------------------------*/
	}

	if (dim == 3) {
		MMSP::grid<3,unsigned long>* grid3=NULL;
		timer = generate<3>(grid3,seeds,nthreads);
		assert(grid3!=NULL);
		#ifdef BGQ
		output_bgq(*grid3, filename);
		#else
		output(*grid3, filename);
		#endif
		#ifndef SILENT
		if (rank==0) std::cout<<"Wrote initial file to "<<filename<<"."<<std::endl;
		#endif
/*------------------Initial tmc----------------------*/
    for(int codx=x0(*grid3, 0); codx < x1(*grid3, 0); codx++) 
      for(int cody=x0(*grid3, 1); cody < x1(*grid3, 1); cody++) 
        for(int codz=x0(*grid3, 2); codz < x1(*grid3, 2); codz++){
          coords[0] = codx;
          coords[1] = cody;
          coords[2] = codz;
          (*grid3).AccessToTmc(coords) = tmc_initial;
          (*grid3).AccessToTmp(coords) = 1.0e6;
        }
/*---------------------------------------------------*/
	}

	return timer;
}

template <int dim>
unsigned long generate(MMSP::grid<dim,unsigned long>*& grid, const char* filename)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	#ifdef MPI_VERSION
	int np = MPI::COMM_WORLD.Get_size();
	#endif

	unsigned long timer=0;
/*------------------Initial tmc----------------------*/
	double tmc_initial = 0;//pow((L_initial/lambda-1)/K1,1.0/n1);
  vector<int> coords (dim,0);
/*---------------------------------------------------*/

	if (dim == 2) {
		grid = new MMSP::grid<dim,unsigned long>(0, 0, dim_x, 0, dim_y);
		(*grid).input(filename, 1, false);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
/*------------------Initial tmc----------------------*/
    for(int codx=x0(*grid, 0); codx < x1(*grid, 0); codx++) 
      for(int cody=x0(*grid, 1); cody < x1(*grid, 1); cody++){
        coords[0] = codx;
        coords[1] = cody;
        (*grid).AccessToTmc(coords) = tmc_initial;
      }
/*---------------------------------------------------*/
	} else if (dim == 3) {
		grid = new MMSP::grid<dim,unsigned long>(0, 0, dim_x, 0, dim_y, 0, dim_z);
		(*grid).input(filename, 1, false);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
/*------------------Initial tmc----------------------*/
    for(int codx=x0(*grid, 0); codx < x1(*grid, 0); codx++) 
      for(int cody=x0(*grid, 1); cody < x1(*grid, 1); cody++) 
        for(int codz=x0(*grid, 2); codz < x1(*grid, 2); codz++){
          coords[0] = codx;
          coords[1] = cody;
          coords[2] = codz;
          (*grid).AccessToTmc(coords) = tmc_initial;
        }
/*---------------------------------------------------*/
	}
	 // double temp[2] = {673,723};
	 // UpdateLocalTmp((*grid), 0.0, &temp[0]);  
	return timer;
}

template <int dim>
unsigned long growthexperiment(MMSP::grid<dim,unsigned long >*& grid, const char* filename)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	#ifdef MPI_VERSION
	int np = MPI::COMM_WORLD.Get_size();
	#endif
	unsigned long timer=0;
	if (dim == 2) {
		grid = new MMSP::grid<dim,unsigned long>(0, 0, dim_x, 0, dim_y);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	} else if (dim == 3) {
		std::cerr << " Initial experiment grain structure input only available for 2D simulation" << std::endl;
	  exit(-1);
	}

/*---------------read grain ID------------------*/
	std::ifstream idtxt(filename, std::ios::in);
	if (!idtxt) {
		std::cerr << "File idtxt error: could not open ";
		std::cerr << filename << "." << std::endl;
	  exit(-1);
  }

  int* ids = new int[dim_x*dim_y];
  int grainID;
  int count = 0;
  int fields = 1;
  while(idtxt>>grainID){
    ids[count] = grainID;
    ++count; 
    fields = (grainID>fields)?grainID:fields;
  }
	MPI::COMM_WORLD.Barrier();
  idtxt.close();
/*----------------------------------------------*/

/*------------------Initial tmc----------------------*/
  double tmc_initial = 0;//pow((L_initial/lambda-1)/K1,1.0/n1);
  vector<int> coords (dim,0);
  if(dim==2){
    for(int codx=x0(*grid, 0); codx < x1(*grid, 0); codx++) 
      for(int cody=x0(*grid, 1); cody < x1(*grid, 1); cody++){
        coords[0] = codx;
        coords[1] = cody;
        (*grid).AccessToTmc(coords) = tmc_initial;
        (*grid).AccessToTmp(coords) = 1.0e6; //set the initial temp to be large enough such that initial physical time = 0.
        (*grid)(coords) = ids[dim_y*codx+cody];
      }
  }
  else if(dim==3){
    std::cerr<<"readin experimental initial grains only support 2D"<<std::endl;
    exit(1);
  }
/*--------------------------------------------------*/
MPI::COMM_WORLD.Barrier();
  delete [] ids;
  ids = NULL;

	return timer;
}

int LargeNearestInteger(int a, int b){
  if(a%b==0) return a/b;
  else return a/b+1;
}

template <int dim> struct flip_index {
	MMSP::grid<dim, unsigned long>* grid;
  int num_of_cells_in_thread;
	int sublattice;
  int num_of_points_to_flip;
  int cell_coord[dim];
  int lattice_cells_each_dimension[dim];
  double Pdenominator;
};

template <int dim> void* flip_index_helper_uniformly( void* s )
{
  srand(time(NULL)); /* seed random number generator */
	flip_index<dim>* ss = static_cast<flip_index<dim>*>(s);
	//  double kT=0.0;
	vector<int> x (dim,0);
  int first_cell_start_coordinates[dim];
  for(int kk=0; kk<dim; kk++) first_cell_start_coordinates[kk] = x0(*(ss->grid), kk);
  for(int i=0; i<dim; i++){
    if(x0(*(ss->grid), i)%2!=0) first_cell_start_coordinates[i]--;
  }
  int cell_coords_selected[dim];
	for (int hh=0; hh<ss->num_of_points_to_flip; hh++) {
	  // choose a random cell to flip
    int cell_numbering_in_thread = rand()%(ss->num_of_cells_in_thread); //choose a cell to flip, from 0 to num_of_cells_in_thread-1
    if(dim==2){
      cell_coords_selected[dim-1]=((ss->cell_coord)[dim-1]+cell_numbering_in_thread)%(ss->lattice_cells_each_dimension)[dim-1];//1-indexed
      cell_coords_selected[0]=(ss->cell_coord)[0]+(((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1]);
    }else if(dim==3){
      cell_coords_selected[dim-1]=((ss->cell_coord)[dim-1]+cell_numbering_in_thread)%(ss->lattice_cells_each_dimension)[dim-1];//1-indexed
      cell_coords_selected[1]=(  (ss->cell_coord)[1]+ ((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1]  )%(ss->lattice_cells_each_dimension)[1];
      cell_coords_selected[0]=(ss->cell_coord)[0]+ ( (ss->cell_coord)[1] + ((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1] ) /(ss->lattice_cells_each_dimension)[1];
    }
    for(int i=0; i<dim; i++){
      x[i]=first_cell_start_coordinates[i]+2*cell_coords_selected[i];
    }
    if(dim==2){
      switch(ss->sublattice){
        case 0:break;// 0,0
        case 1:x[1]++; break; //0,1
        case 2:x[0]++; break; //1,0
        case 3:x[0]++; x[1]++; break; //1,1
      }
    }else if(dim==3){
      switch(ss->sublattice){
        case 0:break;// 0,0,0
        case 1:x[2]++; break; //0,0,1
        case 2:x[1]++; break; //0,1,0
        case 3:x[2]++; x[1]++; break; //0,1,1
        case 4:x[0]++; break; //1,0,0
        case 5:x[2]++; x[0]++; break; //1,0,1
        case 6:x[1]++; x[0]++; break; //1,1,0
        case 7:x[2]++; x[1]++; x[0]++; //1,1,1
      }
    }
    
    bool site_out_of_domain = false;
    for(int i=0; i<dim; i++){
      if(x[i]<x0(*(ss->grid), i) || x[i]>x1(*(ss->grid), i)){
//      if(x[i]<x0(*(ss->grid), i) || x[i]>x1(*(ss->grid), i)-1){
        site_out_of_domain = true;
        break;//break from the for int i loop
      }
    }
    if(site_out_of_domain == true){
      hh--;
      continue; //continue the int hh loop
    }

    int rank = MPI::COMM_WORLD.Get_rank();
    double temperature=273.0;

		unsigned long spin1 = (*(ss->grid))(x);
		// determine neighboring spins
    vector<int> r(dim,0);
    std::vector<unsigned long> neighbors;
    neighbors.clear();
    int number_of_same_neighours = 0;
    if(dim==2){
      for (int i=-1; i<=1; i++) {
        for (int j=-1; j<=1; j++) {
          if(!(i==0 && j==0)){
            r[0] = x[0] + i;
            r[1] = x[1] + j;
            unsigned long spin = (*(ss->grid))(r);

            neighbors.push_back(spin);
            if(spin==spin1) 
              number_of_same_neighours++;
          }
        }
      }
    }else if(dim==3){
		  for (int i=-1; i<=1; i++){
			  for (int j=-1; j<=1; j++){
			    for (int k=-1; k<=1; k++){
            if(!(i==0 && j==0 && k==0)){
				      r[0] = x[0] + i;
				      r[1] = x[1] + j;
				      r[2] = x[2] + k;
				      unsigned long spin = (*(ss->grid))(r);
              neighbors.push_back(spin);
              if(spin==spin1) 
                number_of_same_neighours++;
            }
			    }
        }
		  }
    }
    
    //check if inside a grain
    if(number_of_same_neighours==neighbors.size()){//inside a grain
      continue;//continue for
    }
    //choose a random neighbor spin
    unsigned long spin2 = neighbors[rand()%neighbors.size()];

		if (spin1!=spin2){
			// compute energy change
			double dE = 0.0;
      if(dim==2){
			  for (int i=-1; i<=1; i++){
				  for (int j=-1; j<=1; j++){
            if(!(i==0 && j==0)){
  					  r[0] = x[0] + i;
	  				  r[1] = x[1] + j;
    				  unsigned long spin = (*(ss->grid))(r);
				      dE += 1.0/2*((spin!=spin2)-(spin!=spin1));
            }// if(!(i==0 && j==0))
				  }
        }
      }
      if(dim==3){
			  for (int i=-1; i<=1; i++){ 
				  for (int j=-1; j<=1; j++){ 
    	      for (int k=-1; k<=1; k++){ 
              if(!(i==0 && j==0 && k==0)){
					      r[0] = x[0] + i;
					      r[1] = x[1] + j;
					      r[2] = x[2] + k;
					      unsigned long spin = (*(ss->grid))(r);
					      dE += (spin!=spin2)-(spin!=spin1);
              }
				    }
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

template <int dim> void* flip_index_helper( void* s )
{
  double cos_tolerance=0.999; //cos^-1(0.999)=2.563 deg   cos^-1(0.9999)=0.8103 deg
  int rank = MPI::COMM_WORLD.Get_rank();
  srand((unsigned)time(NULL) + (unsigned)rank); /* seed random number generator */
	flip_index<dim>* ss = static_cast<flip_index<dim>*>(s);
  double kT=0.0;
	vector<int> x (dim,0);
  int first_cell_start_coordinates[dim];
  for(int kk=0; kk<dim; kk++) first_cell_start_coordinates[kk] = x0(*(ss->grid), kk);
  for(int i=0; i<dim; i++){
    if(x0(*(ss->grid), i)%2!=0) first_cell_start_coordinates[i]--;
  }
  int cell_coords_selected[dim];
	for (int hh=0; hh<ss->num_of_points_to_flip; hh++) {
	  // choose a random cell to flip
    int cell_numbering_in_thread = rand()%(ss->num_of_cells_in_thread); //choose a cell to flip, from 0 to num_of_cells_in_thread-1
    if(dim==2){
      cell_coords_selected[dim-1]=((ss->cell_coord)[dim-1]+cell_numbering_in_thread)%(ss->lattice_cells_each_dimension)[dim-1];//1-indexed
      cell_coords_selected[0]=(ss->cell_coord)[0]+(((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1]);
    }else if(dim==3){
      cell_coords_selected[dim-1]=((ss->cell_coord)[dim-1]+cell_numbering_in_thread)%(ss->lattice_cells_each_dimension)[dim-1];//1-indexed
      cell_coords_selected[1]=(  (ss->cell_coord)[1]+ ((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1]  )%(ss->lattice_cells_each_dimension)[1];
      cell_coords_selected[0]=(ss->cell_coord)[0]+ ( (ss->cell_coord)[1] + ((ss->cell_coord)[dim-1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[dim-1] ) /(ss->lattice_cells_each_dimension)[1];
    }
    for(int i=0; i<dim; i++){
      x[i]=first_cell_start_coordinates[i]+2*cell_coords_selected[i];
    }
    if(dim==2){
      switch(ss->sublattice){
        case 0:break;// 0,0
        case 1:x[1]++; break; //0,1
        case 2:x[0]++; break; //1,0
        case 3:x[0]++; x[1]++; break; //1,1
      }
    }else if(dim==3){
      switch(ss->sublattice){
        case 0:break;// 0,0,0
        case 1:x[2]++; break; //0,0,1
        case 2:x[1]++; break; //0,1,0
        case 3:x[2]++; x[1]++; break; //0,1,1
        case 4:x[0]++; break; //1,0,0
        case 5:x[2]++; x[0]++; break; //1,0,1
        case 6:x[1]++; x[0]++; break; //1,1,0
        case 7:x[2]++; x[1]++; x[0]++; //1,1,1
      }
    }

    bool site_out_of_domain = false;
    if(OutsideDomainCheck<dim>(*(ss->grid), &x)){
      hh--;
      continue; //continue the int hh loop
    }

    int rank = MPI::COMM_WORLD.Get_rank();

    double temperature=(*(ss->grid)).AccessToTmp(x);
    double t_mcs = (*(ss->grid)).AccessToTmc(x);
    double Pnumerator = exp(-Q/R/temperature)/pow( pow(L_initial/lambda, m) + K1 * t_mcs, (double)n/m - 1.0) ; //exp(-Q/R/temperature)/pow(t_mcs,n1-1)/pow((1+K1*pow(t_mcs,n1)),n-1);
    double site_selection_probability = Pnumerator/ss->Pdenominator;
	  double rd = double(rand())/double(RAND_MAX);
    if(rd>site_selection_probability){
//      hh--;// no need to guarantee that N times selection is performed in a configurational MC step, so hh is ony need to be the same as in uniform temp case for the highest temp site
      continue;//this site wont be selected
    }

		unsigned long spin1 = (*(ss->grid))(x);
		// determine neighboring spins
    vector<int> r(dim,0);
    std::vector<unsigned long> neighbors;
    neighbors.clear();
    int number_of_same_neighours = 0;
    if(dim==2){
      for (int i=-1; i<=1; i++) {
        for (int j=-1; j<=1; j++) {
          if(!(i==0 && j==0)){
            r[0] = x[0] + i;
            r[1] = x[1] + j;
            if(r[0]<g0(*(ss->grid), 0) || r[0]>=g1(*(ss->grid), 0) || r[1]<g0(*(ss->grid), 1) || r[1]>=g1(*(ss->grid), 1) )
              continue;// neighbour outside the global boundary, skip it. 

            unsigned long spin = (*(ss->grid))(r);
            neighbors.push_back(spin);
            if(spin==spin1) 
              number_of_same_neighours++;
          }
        }
      }
    }else if(dim==3){
		  for (int i=-1; i<=1; i++){
			  for (int j=-1; j<=1; j++){
			    for (int k=-1; k<=1; k++) {
            if(!(i==0 && j==0 && k==0)){
				      r[0] = x[0] + i;
				      r[1] = x[1] + j;
				      r[2] = x[2] + k;

/*----------cancel periodic b.c.-----------*/
	            if(r[0]<g0(*(ss->grid), 0) || r[0]>=g1(*(ss->grid), 0) || r[1]<g0(*(ss->grid), 1) || r[1]>=g1(*(ss->grid), 1) 
                 || r[2]<g0(*(ss->grid), 2) || r[2]>=g1(*(ss->grid), 2))
                continue;// neighbour outside the global boundary, skip it. 
/*----------cancel periodic b.c.-----------*/

				      unsigned long spin = (*(ss->grid))(r);
              neighbors.push_back(spin);
              if(spin==spin1) 
                number_of_same_neighours++;
            }
			    }
        }
		  }
    }
    
    //check if inside a grain
    if(number_of_same_neighours==neighbors.size()){//inside a grain
      continue;//continue for
    }
    //     choose a random neighbor spin
    unsigned long spin2 = neighbors[rand()%neighbors.size()];
		// choose a random spin from Q states
 //       int spin2 = rand()%200;
		if (spin1!=spin2){
			// compute energy change
			double dE = 0.0;
      if(dim==2){
			  for (int i=-1; i<=1; i++){
				  for (int j=-1; j<=1; j++){
            if(!(i==0 && j==0)){
  					  r[0] = x[0] + i;
	  				  r[1] = x[1] + j;

	            if(r[0]<g0(*(ss->grid), 0) || r[0]>=g1(*(ss->grid), 0) || r[1]<g0(*(ss->grid), 1) || r[1]>=g1(*(ss->grid), 1) )
              continue;// neighbour outside the global boundary, skip it. 
    				  unsigned long spin = (*(ss->grid))(r);
				  dE += 1.0/2*((spin!=spin2)-(spin!=spin1));

            }// if(!(i==0 && j==0))
				  }
        }
      }
      if(dim==3){
			  for (int i=-1; i<=1; i++){ 
				  for (int j=-1; j<=1; j++){ 
    	      for (int k=-1; k<=1; k++){ 
              if(!(i==0 && j==0 && k==0)){
					      r[0] = x[0] + i;
					      r[1] = x[1] + j;
					      r[2] = x[2] + k;

/*----------cancel periodic b.c.-----------*/
	            if(r[0]<g0(*(ss->grid), 0) || r[0]>=g1(*(ss->grid), 0) || r[1]<g0(*(ss->grid), 1) || r[1]>=g1(*(ss->grid), 1) 
                 || r[2]<g0(*(ss->grid), 2) || r[2]>=g1(*(ss->grid), 2))
                  continue;// neighbour outside the global boundary, skip it.
/*----------cancel periodic b.c.-----------*/

					      unsigned long spin = (*(ss->grid))(r);
					      dE += 1.0/2*((spin!=spin2)-(spin!=spin1));
              }
				    }
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

template <int dim> bool OutsideDomainCheck(MMSP::grid<dim, unsigned long>& grid, vector<int>* x){
  bool outside_domain=false;
  for(int i=0; i<dim; i++){
    if((*x)[i]<x0(grid, i) || (*x)[i]>x1(grid, i)-1){
      outside_domain=true;
      break;
    }
  }
  return outside_domain;
}

template <int dim> unsigned long update_uniformly(MMSP::grid<dim, int>& grid, int steps, int nthreads)
{
	#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	unsigned int np=0;
//	#ifdef MPI_VERSION
	rank=MPI::COMM_WORLD.Get_rank();
	np=MPI::COMM_WORLD.Get_size();
	MPI::COMM_WORLD.Barrier();
//	#endif

	unsigned long update_timer = 0;

	pthread_t* p_threads = new pthread_t[nthreads];
	flip_index<dim>* mat_para = new flip_index<dim> [nthreads];
	pthread_attr_t attr;
	pthread_attr_init (&attr);

	#ifndef SILENT
	static int iterations = 1;
	if (rank==0) print_progress(0, steps, iterations);
	#endif

/*-----------------------------------------------*/
/*---------------generate cells------------------*/
/*-----------------------------------------------*/
  int dimension_length=0, number_of_lattice_cells=1;
  int lattice_cells_each_dimension[dim];
  for(int i=0; i<dim; i++){
    dimension_length = x1(grid, i)-x0(grid, i);
//    dimension_length = x1(grid, i)-1-x0(grid, i);
    if(x0(grid, 0)%2==0)
      lattice_cells_each_dimension[i] = dimension_length/2+1;
    else
      lattice_cells_each_dimension[i] = 1+LargeNearestInteger(dimension_length,2);
    number_of_lattice_cells *= lattice_cells_each_dimension[i];
  }
  #ifndef SILENT
  if(rank==0){ 
    std::cout<<"lattice_cells_each_dimension are:  ";
    for(int i=0; i<dim; i++)
      std::cout<<lattice_cells_each_dimension[i]<<"  ";
    std::cout<<"\n";
  }
  #endif
//----------assign cells for each pthreads
  int num_of_cells_in_thread = number_of_lattice_cells/nthreads;
  //check if num of the pthread is too large, if so, reduce it.                                                                                                                                         
  if (num_of_cells_in_thread<1) {
    std::cerr<<"ERROR: number of pthread is too large, please reduce it to a value <= "<<number_of_lattice_cells<<std::endl;
    exit(0);
  }

	vector<int> x (dim,0);
	vector<int> x_prim (dim,0);
  int coordinates_of_cell[dim];
  int initial_coordinates[dim];
  int **cell_coord = new int*[nthreads];//record the start coordinates of each pthread domain.
  for(int i=0; i<nthreads; i++){
    cell_coord[i] = new int[dim];
    for(int j=0; j<dim; j++){
      cell_coord[i][j]=0;
    }
  }

  int **num_of_grids_to_flip = new int*[nthreads];
  for(int i=0; i<nthreads; i++){
    num_of_grids_to_flip[i] = new int[( static_cast<int>(pow(2,dim)) )];
    for(int j=0; j<pow(2,dim); j++){
      num_of_grids_to_flip[i][j]=0;
    }
  }

  for(int k=0; k<dim; k++) 
    initial_coordinates[k] = x0(grid, k);
  for(int i=0; i<dim; i++){
    if(x0(grid, i)%2!=0) 
      initial_coordinates[i]--;
  }

  for (int i=0; i<nthreads; i++) {
        int cell_numbering = num_of_cells_in_thread*i; //0-indexed, celling_numbering is the start cell numbering
        if(dim==2){
          cell_coord[i][dim-1]=cell_numbering%lattice_cells_each_dimension[dim-1];//0-indexed
          cell_coord[i][0]=(cell_numbering/lattice_cells_each_dimension[dim-1]);
	        if(cell_coord[i][0]>=lattice_cells_each_dimension[0]){
	          std::cerr<<"the cell coordinates is wrong!"<<std::endl;
	          exit(1);
	        }
        }else if(dim==3){
          cell_coord[i][dim-1]=cell_numbering%lattice_cells_each_dimension[dim-1];//0-indexed
          cell_coord[i][1]=(cell_numbering/lattice_cells_each_dimension[dim-1])%lattice_cells_each_dimension[1];
          cell_coord[i][0]=(cell_numbering/lattice_cells_each_dimension[dim-1])/lattice_cells_each_dimension[1];
	        if(cell_coord[i][0]>=lattice_cells_each_dimension[0]){
	          std::cerr<<"the cell coordinates is wrong!"<<std::endl;
	          exit(1);
	        }
        }

    mat_para[i].grid = &grid;
    if(i==(nthreads-1)) 
      mat_para[i].num_of_cells_in_thread = number_of_lattice_cells - num_of_cells_in_thread*(nthreads-1);
    else 
      mat_para[i].num_of_cells_in_thread = num_of_cells_in_thread;

    #ifndef SILENT
    if(rank==0) std::cout<<"num_of_cells_in_thread is "<<mat_para[i].num_of_cells_in_thread<<" in thread "<<i<<"\n";
    #endif

    for(int k=0; k<dim; k++) 
      mat_para[i].lattice_cells_each_dimension[k]=lattice_cells_each_dimension[k];

//    mat_para[i].temperature_along_x = temperature_along_x;
//    mat_para[i].grain_orientations = grain_orientations;

    for(int j=0; j<mat_para[i].num_of_cells_in_thread; j++){
      int start_cell_numbering = num_of_cells_in_thread*i;
      if(dim==2){
        coordinates_of_cell[dim-1]=(start_cell_numbering+j)%lattice_cells_each_dimension[dim-1];//0-indexed
        coordinates_of_cell[0]=(start_cell_numbering+j)/lattice_cells_each_dimension[dim-1];
      }else if(dim==3){
        coordinates_of_cell[dim-1]=(start_cell_numbering+j)%lattice_cells_each_dimension[dim-1];//0-indexed
        coordinates_of_cell[1]=((start_cell_numbering+j)/lattice_cells_each_dimension[dim-1])%lattice_cells_each_dimension[1];
        coordinates_of_cell[0]=((start_cell_numbering+j)/lattice_cells_each_dimension[dim-1])/lattice_cells_each_dimension[1];
      }
      for(int ii=0; ii<dim; ii++){
        x[ii]=initial_coordinates[ii]+2*coordinates_of_cell[ii];
      }

      if(dim==2){
        x_prim = x;
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][0]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1; //0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][1]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1; //1,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][2]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1;
        x_prim[1]=x[1]+1; //1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][3]+=1;
      }else if(dim==3){
        x_prim = x;//0,0,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][0]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1; //0,0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][1]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1; //0,1,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][2]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[1]=x[1]+1; //0,1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][3]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1; //1,0,0 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][4]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[0]=x[0]+1; //1,0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][5]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1;
        x_prim[0]=x[0]+1; //1,1,0 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][6]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[1]=x[1]+1;
        x_prim[0]=x[0]+1; //1,1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][7]+=1;
      }
    }// for int j 
  }//for int i

	for (int step=0; step<steps; step++){
		unsigned long start = rdtsc();
    int num_of_sublattices=0;
    if(dim==2) num_of_sublattices = 4; 
    else if(dim==3) num_of_sublattices = 8;
		for (int sublattice=0; sublattice < num_of_sublattices; sublattice++) {
			for (int i=0; i!= nthreads ; i++) {
				mat_para[i].sublattice=sublattice;
				mat_para[i].num_of_points_to_flip=num_of_grids_to_flip[i][sublattice];
        for(int k=0; k<dim; k++) mat_para[i].cell_coord[k]=cell_coord[i][k];
				pthread_create(&p_threads[i], &attr, flip_index_helper_uniformly<dim>, (void*) &mat_para[i] );
			}//loop over threads

			for (int ii=0; ii!= nthreads ; ii++)
				pthread_join(p_threads[ii], NULL);

			#ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
			#endif

//			ghostswap(grid, sublattice); // once looped over a "color", ghostswap.
			ghostswap(grid); // once looped over a "color", ghostswap.
      #ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
                #endif

		}//loop over color
		#ifndef SILENT
		if (rank==0) print_progress(step+1, steps, iterations);
		#endif
		update_timer += rdtsc()-start;
	}//loop over step
	#ifndef SILENT
	++iterations;
	#endif

  for(int i=0; i<nthreads; i++){
    delete [] num_of_grids_to_flip[i];
    num_of_grids_to_flip[i]=NULL;
    delete [] cell_coord[i];
    cell_coord[i]=NULL;
  }
  delete [] num_of_grids_to_flip; 
  num_of_grids_to_flip=NULL; 
  delete [] cell_coord;
  cell_coord=NULL;
//	delete [] grain_orientations;
//	grain_orientations=NULL;
//	delete [] temperature_along_x;
//	temperature_along_x=NULL;
	delete [] p_threads;
	p_threads=NULL;
	delete [] mat_para;
	mat_para=NULL;
	unsigned long total_update_time=update_timer;
	#ifdef MPI_VERSION
	MPI::COMM_WORLD.Allreduce(&update_timer, &total_update_time, 1, MPI_UNSIGNED_LONG, MPI_SUM);
	#endif
	return total_update_time/np; // average update time
}

template <int dim> double PdenominatorMax(MMSP::grid<dim, unsigned long>& grid, double& tmc_at_PdenominatorMax, double& tmp_at_PdenominatorMax){
   double Pdenominator_max = 0.0;
   vector<int> coords (dim,0);
   if(dim==2){
       for(int codx=x0(grid, 0); codx < x1(grid, 0); codx++) 
         for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
           coords[0] = codx;
           coords[1] = cody;
           double temperature = grid.AccessToTmp(coords);
           double tmc_local = grid.AccessToTmc(coords);
           double Pdenominator = exp(-Q/R/temperature)/pow( pow(L_initial/lambda, m) + K1 * tmc_local, (double)n/m - 1.0); //exp(-Q/R/temperature)/pow(tmc_local,n1-1)/pow((1+K1*pow(tmc_local,n1)),n-1);
           if(Pdenominator > Pdenominator_max){
             Pdenominator_max = Pdenominator;
             tmc_at_PdenominatorMax = grid.AccessToTmc(coords);
             tmp_at_PdenominatorMax = temperature;
           }
         }
   }
   else if(dim==3){
     for(int codx=x0(grid, 0); codx < x1(grid, 0); codx++) 
       for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++) 
         for(int codz=x0(grid, 2); codz < x1(grid, 2); codz++){
           coords[0] = codx;
           coords[1] = cody;
           coords[2] = codz;
           double temperature = grid.AccessToTmp(coords);
           double tmc_local = grid.AccessToTmc(coords);
           double Pdenominator = exp(-Q/R/temperature)/pow( pow(L_initial/lambda, m) + K1 * tmc_local, (double)n/m - 1.0); //exp(-Q/R/temperature)/pow(tmc_local,n1-1)/pow((1+K1*pow(tmc_local,n1)),n-1);                                                                   
           if(Pdenominator > Pdenominator_max){
             Pdenominator_max = Pdenominator;
             tmc_at_PdenominatorMax = grid.AccessToTmc(coords);
             tmp_at_PdenominatorMax = temperature;
           }
        }
   }
   return Pdenominator_max;
}

/*template <int dim> void UpdateLocalTmc(MMSP::grid<dim, unsigned long>& grid, double t_inc){
   vector<int> coords (dim,0);
   if(dim==2){
       for(int codx=x0(grid, 0); codx < x1(grid, 0); codx++) 
         for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
           coords[0] = codx;
           coords[1] = cody;
           double tmc_local = grid.AccessToTmc(coords);
           long double exp_eqn_rhs = pow( pow(lambda, m) * K1 * tmc_local + pow(L_initial, m), (double)n/m) - pow(L0, n); //pow(lambda*(1+K1*pow(tmc_local, n1)), n) - pow(L0, n);
           long double temperature = grid.AccessToTmp(coords);
           exp_eqn_rhs += K_*t_inc*exp(-Q/R/temperature);
           grid.AccessToTmc(coords) = (pow(exp_eqn_rhs + pow(L0, n), (double)m/n) - pow(L_initial, m)) / pow(lambda, m) / K1;//pow(1.0/K1/lambda*pow(exp_eqn_rhs+pow(L0, n),1.0/n)-1.0/K1, 1.0/n1);  
         }
   }
   else if(dim==3){
     for(int codx=x0(grid, 0); codx < x1(grid, 0); codx++)  
         for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++) 
           for(int codz=x0(grid, 2); codz < x1(grid, 2); codz++){
             coords[0] = codx;
             coords[1] = cody;
             coords[2] = codz;
             double tmc_local = grid.AccessToTmc(coords);
             long double exp_eqn_rhs = pow( pow(lambda, m) * K1 * tmc_local + pow(L_initial, m), (double)n/m) - pow(L0, n); //pow(lambda*(1+K1*pow(tmc_local, n1)), n) - pow(L0, n);
             long double temperature = grid.AccessToTmp(coords);
             exp_eqn_rhs += K_*t_inc*exp(-Q/R/temperature);
             grid.AccessToTmc(coords) = (pow(exp_eqn_rhs + pow(L0, n), (double)m/n) - pow(L_initial, m)) / pow(lambda, m) / K1;//pow(1.0/K1/lambda*pow(exp_eqn_rhs+pow(L0, n),1.0/n)-1.0/K1, 1.0/n1);  
           }
   }
}*/

/*template <int dim> void calCulateGrainSize(MMSP::grid<dim, unsigned long>& grid, unsigned long &number_of_grains){
   vector<int> coords (dim,0);
   if(dim==2){
      int midx = (g0(grid, 0)+g1(grid, 0))/2;
      int midy = (g0(grid, 1)+g1(grid, 1))/2;
      if(midx>=x0(grid, 0) && midx<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = midx;
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.15*midx>=x0(grid, 0) && 0.15*midx<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.15*midx;
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.35*midx>=x0(grid, 0) && 0.35*midx<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.35*midx;
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.65*midx>=x0(grid, 0) && 0.65*midx<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.65*midx;
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.85*midx>=x0(grid, 0) && 0.85*midx<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.85*midx;
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
   }
   else if(dim==3){
      int midx = (g0(grid, 0)+g1(grid, 0))/2;
      int midy = (g0(grid, 1)+g1(grid, 1))/2;
      int midz = (g0(grid, 2)+g1(grid, 2))/2;
      if(midx>=x0(grid, 0) && midx<x1(grid, 0) && midz>=x0(grid, 2) && midz<x1(grid, 2)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = midx;
          coords[1] = cody;
          coords[2] = midz;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.25*midx>=x0(grid, 0) && 0.25*midx<x1(grid, 0) && 0.25*midz>=x0(grid, 2) && 0.25*midz<x1(grid, 2)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.25*midx;
          coords[1] = cody;
          coords[2] = 0.25*midz;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.75*midx>=x0(grid, 0) && 0.75*midx<x1(grid, 0) && 0.75*midz>=x0(grid, 2) && 0.75*midz<x1(grid, 2)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.75*midx;
          coords[1] = cody;
          coords[2] = 0.75*midz;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.75*midx>=x0(grid, 0) && 0.75*midx<x1(grid, 0) && 0.25*midz>=x0(grid, 2) && 0.25*midz<x1(grid, 2)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.75*midx;
          coords[1] = cody;
          coords[2] = 0.25*midz;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.25*midx>=x0(grid, 0) && 0.25*midx<x1(grid, 0) && 0.75*midz>=x0(grid, 2) && 0.75*midz<x1(grid, 2)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.25*midx;
          coords[1] = cody;
          coords[2] = 0.75*midz;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
   }
}*/

/*template <int dim> void calCulateGrainSize1(MMSP::grid<dim, unsigned long>& grid, unsigned long &number_of_grains){
   vector<int> coords (dim,0);
   if(dim==2){
      if(0.08*g1(grid, 0)>=x0(grid, 0) && 0.08*g1(grid, 0)<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.08*g1(grid, 0);
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.16*g1(grid, 0)>=x0(grid, 0) && 0.16*g1(grid, 0)<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.16*g1(grid, 0);
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.24*g1(grid, 0)>=x0(grid, 0) && 0.24*g1(grid, 0)<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.24*g1(grid, 0);
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.32*g1(grid, 0)>=x0(grid, 0) && 0.32*g1(grid, 0)<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.32*g1(grid, 0);
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.40*g1(grid, 0)>=x0(grid, 0) && 0.40*g1(grid, 0)<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.40*g1(grid, 0);
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
   }
}
*/
/*template <int dim> void calCulateGrainSize2(MMSP::grid<dim, unsigned long>& grid, unsigned long &number_of_grains){
   vector<int> coords (dim,0);
   if(dim==2){
      if(0.60*g1(grid, 0)>=x0(grid, 0) && 0.60*g1(grid, 0)<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.60*g1(grid, 0);
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.68*g1(grid, 0)>=x0(grid, 0) && 0.68*g1(grid, 0)<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.68*g1(grid, 0);
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.76*g1(grid, 0)>=x0(grid, 0) && 0.76*g1(grid, 0)<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.76*g1(grid, 0);
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.84*g1(grid, 0)>=x0(grid, 0) && 0.84*g1(grid, 0)<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.84*g1(grid, 0);
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
      if(0.92*g1(grid, 0)>=x0(grid, 0) && 0.92*g1(grid, 0)<x1(grid, 0)){
        unsigned long index_record = -1;
        for(int cody=x0(grid, 1); cody < x1(grid, 1); cody++){
          coords[0] = 0.92*g1(grid, 0);
          coords[1] = cody;
          if(grid(coords)!=index_record){
            index_record = grid(coords);
            number_of_grains = number_of_grains + 1;       
          }
        }
        if(x1(grid, 1) < g1(grid, 1)){
          coords[1] = x1(grid, 1);
          if(grid(coords) == index_record && number_of_grains != 0){
            number_of_grains = number_of_grains - 1; 
          }
        }
      }
   }
}*/

template <int dim> unsigned long update(MMSP::grid<dim, unsigned long>& grid, int steps, int steps_finished, int nthreads, int step_to_nonuniform, long double &physical_time, double* temp)
{
	#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	unsigned int np=0;
//	#ifdef MPI_VERSION
	rank=MPI::COMM_WORLD.Get_rank();
	np=MPI::COMM_WORLD.Get_size();
	MPI::COMM_WORLD.Barrier();
//	#endif

	unsigned long update_timer = 0;

	pthread_t* p_threads = new pthread_t[nthreads];
	flip_index<dim>* mat_para = new flip_index<dim> [nthreads];
	pthread_attr_t attr;
	pthread_attr_init (&attr);

	#ifndef SILENT
	static int iterations = 1;
	if (rank==0) print_progress(0, steps, iterations);
	#endif
	ghostswap(grid); 
/*
  int edge = 1000;
  int number_of_fields = static_cast<int>(float(edge*edge)/(M_PI*10.*10.)); // average grain is a disk of radius 10
  #ifdef MPI_VERSION
	while (number_of_fields % np) --number_of_fields;
  #endif
*/
/*
  EulerAngles *grain_orientations = new EulerAngles[200];
  if(rank==0){
    ReadData(grain_orientations);
  }
	MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Bcast(grain_orientations, 200*3, MPI_DOUBLE, 0);
*/
/*-----------------------------------------------*/
/*---------------generate cells------------------*/
/*-----------------------------------------------*/
  int dimension_length=0, number_of_lattice_cells=1;
  int lattice_cells_each_dimension[dim];
  for(int i=0; i<dim; i++){
    dimension_length = x1(grid, i)-x0(grid, i);
//    dimension_length = x1(grid, i)-1-x0(grid, i);
    if(x0(grid, 0)%2==0)
      lattice_cells_each_dimension[i] = dimension_length/2+1;
    else
      lattice_cells_each_dimension[i] = 1+LargeNearestInteger(dimension_length,2);
    number_of_lattice_cells *= lattice_cells_each_dimension[i];
  }
  #ifndef SILENT
  if(rank==0){ 
    std::cout<<"lattice_cells_each_dimension are:  ";
    for(int i=0; i<dim; i++)
      std::cout<<lattice_cells_each_dimension[i]<<"  ";
    std::cout<<"\n";
  }
  #endif
//----------assign cells for each pthreads
  int num_of_cells_in_thread = number_of_lattice_cells/nthreads;
  //check if num of the pthread is too large, if so, reduce it.                                                                                                                                         
  if (num_of_cells_in_thread<1) {
    std::cerr<<"ERROR: number of pthread is too large, please reduce it to a value <= "<<number_of_lattice_cells<<std::endl;
    exit(0);
  }


	vector<int> x (dim,0);
	vector<int> x_prim (dim,0);
  int coordinates_of_cell[dim];
  int initial_coordinates[dim];
  int **cell_coord = new int*[nthreads];//record the start coordinates of each pthread domain.
  for(int i=0; i<nthreads; i++){
    cell_coord[i] = new int[dim];
    for(int j=0; j<dim; j++){
      cell_coord[i][j]=0;
    }
  }

  int **num_of_grids_to_flip = new int*[nthreads];
  for(int i=0; i<nthreads; i++){
    num_of_grids_to_flip[i] = new int[( static_cast<int>(pow(2,dim)) )];
    for(int j=0; j<pow(2,dim); j++){
      num_of_grids_to_flip[i][j]=0;
    }
  }

  for(int k=0; k<dim; k++) 
    initial_coordinates[k] = x0(grid, k);
  for(int i=0; i<dim; i++){
    if(x0(grid, i)%2!=0) 
      initial_coordinates[i]--;
  }

  for (int i=0; i<nthreads; i++) {
        int cell_numbering = num_of_cells_in_thread*i; //0-indexed, celling_numbering is the start cell numbering
        if(dim==2){
          cell_coord[i][dim-1]=cell_numbering%lattice_cells_each_dimension[dim-1];//0-indexed
          cell_coord[i][0]=(cell_numbering/lattice_cells_each_dimension[dim-1]);
	        if(cell_coord[i][0]>=lattice_cells_each_dimension[0]){
	          std::cerr<<"the cell coordinates is wrong!"<<std::endl;
	          exit(1);
	        }
        }else if(dim==3){
          cell_coord[i][dim-1]=cell_numbering%lattice_cells_each_dimension[dim-1];//0-indexed
          cell_coord[i][1]=(cell_numbering/lattice_cells_each_dimension[dim-1])%lattice_cells_each_dimension[1];
          cell_coord[i][0]=(cell_numbering/lattice_cells_each_dimension[dim-1])/lattice_cells_each_dimension[1];
	        if(cell_coord[i][0]>=lattice_cells_each_dimension[0]){
	          std::cerr<<"the cell coordinates is wrong!"<<std::endl;
	          exit(1);
	        }
        }

    mat_para[i].grid = &grid;
    if(i==(nthreads-1)) 
      mat_para[i].num_of_cells_in_thread = number_of_lattice_cells - num_of_cells_in_thread*(nthreads-1);
    else 
      mat_para[i].num_of_cells_in_thread = num_of_cells_in_thread;

    #ifndef SILENT
    if(rank==0) std::cout<<"num_of_cells_in_thread is "<<mat_para[i].num_of_cells_in_thread<<" in thread "<<i<<"\n";
    #endif

    for(int k=0; k<dim; k++) 
      mat_para[i].lattice_cells_each_dimension[k]=lattice_cells_each_dimension[k];

    for(int j=0; j<mat_para[i].num_of_cells_in_thread; j++){
      int start_cell_numbering = num_of_cells_in_thread*i;
      if(dim==2){
        coordinates_of_cell[dim-1]=(start_cell_numbering+j)%lattice_cells_each_dimension[dim-1];//0-indexed
        coordinates_of_cell[0]=(start_cell_numbering+j)/lattice_cells_each_dimension[dim-1];
      }else if(dim==3){
        coordinates_of_cell[dim-1]=(start_cell_numbering+j)%lattice_cells_each_dimension[dim-1];//0-indexed
        coordinates_of_cell[1]=((start_cell_numbering+j)/lattice_cells_each_dimension[dim-1])%lattice_cells_each_dimension[1];
        coordinates_of_cell[0]=((start_cell_numbering+j)/lattice_cells_each_dimension[dim-1])/lattice_cells_each_dimension[1];
      }
      for(int ii=0; ii<dim; ii++){
        x[ii]=initial_coordinates[ii]+2*coordinates_of_cell[ii];
      }

      if(dim==2){
        x_prim = x;
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][0]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1; //0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][1]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1; //1,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][2]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1;
        x_prim[1]=x[1]+1; //1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][3]+=1;
      }else if(dim==3){
        x_prim = x;//0,0,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][0]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1; //0,0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][1]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1; //0,1,0
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][2]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[1]=x[1]+1; //0,1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][3]+=1;

        x_prim = x;
        x_prim[0]=x[0]+1; //1,0,0 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][4]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[0]=x[0]+1; //1,0,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][5]+=1;

        x_prim = x;
        x_prim[1]=x[1]+1;
        x_prim[0]=x[0]+1; //1,1,0 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][6]+=1;

        x_prim = x;
        x_prim[2]=x[2]+1;
        x_prim[1]=x[1]+1;
        x_prim[0]=x[0]+1; //1,1,1 
        if(!OutsideDomainCheck<dim>(grid, &x_prim)) num_of_grids_to_flip[i][7]+=1;
      }
    }// for int j 
  }//for int i

	for (int step=0; step<steps; step++){

    double tmc_at_PdenominatorMax = 0.0;
    double tmp_at_PdenominatorMax = 0.0;
    double Pdenominator_max_partition = PdenominatorMax(grid, tmc_at_PdenominatorMax, tmp_at_PdenominatorMax); 
    MPI::COMM_WORLD.Barrier();
    double Pdenominator_max_global = 0.0;

    MPI::COMM_WORLD.Allreduce(&Pdenominator_max_partition, &Pdenominator_max_global, 1, MPI_DOUBLE, MPI_MAX);
    MPI::COMM_WORLD.Barrier();

    double tmc_at_PdenominatorMax_global = 0.0;
    double tmp_at_PdenominatorMax_global = 0.0;
    if(Pdenominator_max_partition != Pdenominator_max_global){
      tmc_at_PdenominatorMax = 0.0;
      tmp_at_PdenominatorMax = 0.0;
    }
    MPI::COMM_WORLD.Allreduce(&tmc_at_PdenominatorMax, &tmc_at_PdenominatorMax_global, 1, MPI_DOUBLE, MPI_MAX);
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Allreduce(&tmp_at_PdenominatorMax, &tmp_at_PdenominatorMax_global, 1, MPI_DOUBLE, MPI_MAX);
    MPI::COMM_WORLD.Barrier();

    double t_inc = (pow(pow(lambda, m) * K1 * (tmc_at_PdenominatorMax_global + 1) + pow(L_initial, m), (double)n/m) - pow(pow(lambda, m) * K1 * tmc_at_PdenominatorMax_global + pow(L_initial, m), (double)n/m)) / K_/exp(-Q/R/tmp_at_PdenominatorMax_global); //(pow(lambda*(1+K1*pow(tmc_at_PdenominatorMax_global+1, n1)), n) - 
      // pow(lambda*(1+K1*pow(tmc_at_PdenominatorMax_global, n1)), n) )/K_/exp(-Q/R/tmp_at_PdenominatorMax_global);

		unsigned long start = rdtsc();
    int num_of_sublattices=0;
    if(dim==2) num_of_sublattices = 4; 
    else if(dim==3) num_of_sublattices = 8;
		for (int sublattice=0; sublattice < num_of_sublattices; sublattice++) {
			for (int i=0; i!= nthreads ; i++) {
				mat_para[i].sublattice=sublattice;
				mat_para[i].num_of_points_to_flip=num_of_grids_to_flip[i][sublattice];
        mat_para[i].Pdenominator = Pdenominator_max_global;
        for(int k=0; k<dim; k++) mat_para[i].cell_coord[k]=cell_coord[i][k];
				pthread_create(&p_threads[i], &attr, flip_index_helper<dim>, (void*) &mat_para[i] );
			}//loop over threads

			for (int ii=0; ii!= nthreads ; ii++)
				pthread_join(p_threads[ii], NULL);

			#ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
			#endif

			ghostswap(grid, sublattice); // once looped over a "color", ghostswap.
//			ghostswap(grid); // once looped over a "color", ghostswap.
            #ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
                #endif
		}//loop over color


	  MPI::COMM_WORLD.Barrier();
    physical_time += t_inc;
	  MPI::COMM_WORLD.Barrier();
   //  UpdateLocalTmc(grid, t_inc);
	  // MPI::COMM_WORLD.Barrier();
   // UpdateLocalTmp(grid, physical_time, temp);
	  // MPI::COMM_WORLD.Barrier();

//     unsigned long number_of_grains = 0;
//     calCulateGrainSize1(grid, number_of_grains);
//     unsigned long total_number_of_grains = 0;
//     MPI::COMM_WORLD.Allreduce(&number_of_grains, &total_number_of_grains, 1, MPI_UNSIGNED_LONG, MPI_SUM);
//     MPI::COMM_WORLD.Barrier();

//     double grain_size1 =  L0 + lambda*(5.0 * (dim_y / total_number_of_grains));

//     number_of_grains = 0;
//     calCulateGrainSize2(grid, number_of_grains);
//     total_number_of_grains = 0;
//     MPI::COMM_WORLD.Allreduce(&number_of_grains, &total_number_of_grains, 1, MPI_UNSIGNED_LONG, MPI_SUM);
//     MPI::COMM_WORLD.Barrier();

//     double grain_size2 = L0 + lambda* (5.0 * (dim_y / total_number_of_grains));

// if(rank==0)
// std::cout<< "  "<<tmc_at_PdenominatorMax_global<< "  "<<grain_size1<<"  "<<grain_size2<< std::endl;
// 		update_timer += rdtsc()-start;
	}//loop over step

	#ifndef SILENT
	++iterations;
	#endif

  for(int i=0; i<nthreads; i++){
    delete [] num_of_grids_to_flip[i];
    num_of_grids_to_flip[i]=NULL;
    delete [] cell_coord[i];
    cell_coord[i]=NULL;
  }
  delete [] num_of_grids_to_flip; 
  num_of_grids_to_flip=NULL; 
  delete [] cell_coord;
  cell_coord=NULL;
	delete [] p_threads;
	p_threads=NULL;
	delete [] mat_para;
	mat_para=NULL;

	unsigned long total_update_time=update_timer;
	#ifdef MPI_VERSION
	MPI::COMM_WORLD.Allreduce(&update_timer, &total_update_time, 1, MPI_UNSIGNED_LONG, MPI_SUM);
  MPI::COMM_WORLD.Barrier();
	#endif
	return total_update_time/np; // average update time
}

}


#endif

#include"MMSP.main.hpp"

