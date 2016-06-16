/*=========================================\\
|| Daisy World  Model Class implementation ||
\\=========================================*/
#include "DaisyWorld.h"

// Static class parameters
const double DaisyWorld::sigmaB;
const double DaisyWorld::S0;
const double DaisyWorld::T0;
const double DaisyWorld::DT;
const double DaisyWorld::CT;
const double DaisyWorld::A0;
const double DaisyWorld::gamma;
const double DaisyWorld::Tmin;
const double DaisyWorld::Tmax;
const double DaisyWorld::Tideal;

// Main constructor
// Requires grid size (N), solar luminosity factor (solar_lumin),
// mutation rate (mutation_rate) and noise amplitude (noise_sigma)
DaisyWorld::DaisyWorld(int p_N, double p_solar_lumin, double p_mutation_rate, double p_noise_sigma) {

  int i, j;

  N = p_N;
  solar_lumin = p_solar_lumin;
  mutation_rate = p_mutation_rate;
  noise_sigma = p_noise_sigma;
  printf("noise = %f\n", noise_sigma);

  // Allocate data arrays (and buffers)
  vegetation  = new bool*[N];
  temperature = new double*[N];
  albedo      = new double*[N];
  vegetation_buffer  = new bool*[N];
  temperature_buffer = new double*[N];
  albedo_buffer      = new double*[N];  
  for (i = 0; i < N; i++) {
    vegetation[i]  = new bool[N];
    temperature[i] = new double[N];
    albedo[i]      = new double[N];
    vegetation_buffer[i]  = new bool[N];
    temperature_buffer[i] = new double[N];
    albedo_buffer[i]      = new double[N];
  }

  // Misc stuff
  num_cells = N*N;
  growth_constant = 4.0/(Tmax-Tmin)/(Tmax-Tmin);
  
}

// =========================================

// Returns a random real uniformly distributed in [0,1)
double DaisyWorld::random_real() {
  return ((double)rand())/((double)RAND_MAX);
}

// =========================================

// Reset the simulation, using a newly generated randomized seed
void DaisyWorld::reset() {
  srand(time(NULL));
  global_seed = rand();
  reset(global_seed);
}

// =========================================

// Reset the simulation to its initial state, using the provided seed
void DaisyWorld::reset(int seed) {

  // Seed RNG
  srand(seed);

  // Randomly seed daisies -- vegetation coverage is 100%, local
  // temperatures are set to the thermodynamic equlibrium temperature
  Teq = pow(solar_lumin*S0*(1.0-A0)/sigmaB, 0.25);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      vegetation[i][j] = true;
      albedo[i][j] = random_real();
      temperature[i][j] = Teq;
    }
  }

  global_temperature = compute_global_temperature();
  global_albedo = compute_global_albedo();
  global_vegetation = compute_global_vegetation();
  it = 0;

}

// =========================================

// Compute global average temperature
double DaisyWorld::compute_global_temperature() {
  double temp = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      temp += temperature[i][j];
    }
  }
  return temp / ((double) num_cells);
}

// =========================================

// Compute global average albedo
double DaisyWorld::compute_global_albedo() {
  double alb = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      alb += albedo[i][j];
    }
  }
  return alb / ((double) num_cells);
}

// =========================================

// Compute global vegetation coverage
double DaisyWorld::compute_global_vegetation() {
  double veg = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (vegetation[i][j]) veg += 1.0;
    }
  }
  return veg / ((double) num_cells);
}

// =========================================

// Perform one iteration of the model
void DaisyWorld::time_update() {

  // Apply daisy dynamics
  update_vegetation();

  // Solve temperature dynamics
  update_temperature();

  // Update global temperature, albedo and vegetation coverage
  global_temperature = compute_global_temperature();
  global_albedo = compute_global_albedo();
  global_vegetation = compute_global_vegetation();

  it++;

}

// =========================================

// Applies the cellular-automaton rules for daisy dynamics,
// updating coverage[][] and albedo[][]
void DaisyWorld::update_vegetation() {

  int neigh, ni, nj;
  bool verbose = false;

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {

      if (vegetation[i][j]) {

        // If cell has vegetation, kill it with prob. gamma
        if (verbose) printf("Cell %d,%d has vegetation", i, j);
        if (random_real() <= gamma) {
          if (verbose) printf(" -> died\n");
          vegetation_buffer[i][j] = false;
          albedo_buffer[i][j] = A0; 
        } else {
          if (verbose) printf(" -> survived\n");
          vegetation_buffer[i][j] = true;
          albedo_buffer[i][j] = albedo[i][j]; 
        } 

      } else {

        // If the cell is without vegetation, randomly pick a neighbor
        if (verbose) printf("Cell %d,%d is empty", i, j);
        neigh = rand() % 8;
        get_neigh_coords(i, j, neigh, ni, nj);
        if (verbose) printf("-> neighbor %d,%d ", ni, nj);

        if (vegetation[ni][nj]) {

          // If neighbor has vegetation, growth occurs with prob. beta
          if (verbose) printf("has vegetation ");
          if (verbose) printf("-> (beta=%.3f) ", compute_beta(temperature[ni][nj])); 
          if (random_real() <= compute_beta(temperature[ni][nj])) {
            if (verbose) printf("spawned growth\n");
            vegetation_buffer[i][j] = true;
            albedo_buffer[i][j] = child_albedo(albedo[ni][nj]);
          } else {
            if (verbose) printf("did not span growth\n");
            vegetation_buffer[i][j] = false;
            albedo_buffer[i][j] = A0; 
          }
          
        } else {

          // If neighbor does not have vegetation, cell remains empty
          if (verbose) printf("is also empty\n");
          vegetation_buffer[i][j] = false;
          albedo_buffer[i][j] = A0; 

        }

      }

    }
  }

  // Swap buffer references
  bool** tmp_bool;
  tmp_bool = vegetation;
  vegetation = vegetation_buffer;
  vegetation_buffer = tmp_bool;
  double** tmp_double;
  tmp_double = albedo;
  albedo = albedo_buffer;
  albedo_buffer = tmp_double;


}

// =========================================

// Updates the temperatue distribution
// This is done by solving the energy balance PDE through a simple
// explicit difference scheme, with the albedo distribution held fixed
// The global constant CT controls the timescale of the temperature
// dynamics relative to that of the vegetation dynamics; the larger the
// value, the longer the timescale of the temperature dynamics
void DaisyWorld::update_temperature() {

  int ip, im, jp, jm;

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {

      // Neighbor coordinates
      im = i-1; ip = i+1;
      apply_topology(im, ip);
      jm = j-1; jp = j+1;
      apply_topology(jm, jp);

      // Explicit first-order forward-time centered-space finite
      // difference scheme
      temperature_buffer[i][j] = temperature[i][j] + 1.0/CT * (
      + DT * (temperature[ip][j] + temperature[im][j]
        + temperature[i][jp] + temperature[i][jm]
        - 4*temperature[i][j])
      - sigmaB * pow(temperature[i][j], 4)
      + solar_lumin * S0 * (1 - albedo[i][j])
      ) + noise_term();

    }
  }

  // Swap buffer references
  double** tmp_double;
  tmp_double = temperature;
  temperature = temperature_buffer;
  temperature_buffer = tmp_double;

}

// =========================================

// Returns the albedo of a child daisy patch based on the albedo of
// the progenitor patch
// Possible mutation is determined here, based on the mutation rate (=R)
double DaisyWorld::child_albedo(double prog_albedo) {
  if (mutation_rate == 0) {
    return prog_albedo;
  } else {
    // Mutation model: Uniform albedo drift in [-R,R]
    return prog_albedo + mutation_rate*(2*random_real() - 1.0);
  }
}

// =========================================

// Computes the growth probability beta, as a function of temperature
// This is computed as:  
//   growth_constant * (T-Tmin) * (Tmax-T)
// where
//   growth_constant = 4.0/(Tmax-Tmin)^2
double DaisyWorld::compute_beta(double T) {
  if ((T < Tmin) || (T > Tmax)) {
    return 0;
  } else {
    return growth_constant*(T-Tmin)*(Tmax-T);
  }
}

// =========================================

// Returns the coordinates of the neighbour of cell (i,j) numbered
// neigh using the provided ni and nj variable references
// If i and j are thought to increase to the right and up, respectively,
// then neighbor numbering is as follows:
//  0  1  2
//  7  .  3
//  6  5  4
void DaisyWorld::get_neigh_coords(int i, int j, int neigh, int& ni, int& nj) {

  switch (neigh) {
    case 0: ni = i-1; nj = j+1; break;
    case 1: ni = i  ; nj = j+1; break;
    case 2: ni = i+1; nj = j+1; break;
    case 3: ni = i+1; nj = j  ; break;
    case 4: ni = i+1; nj = j-1; break;
    case 5: ni = i  ; nj = j-1; break;
    case 6: ni = i-1; nj = j-1; break;
    case 7: ni = i-1; nj = j  ; break;
    default: printf("WRONG NEIGHBOR NUMBER: %d!\n", neigh); exit(1);
  }
  
  apply_topology(ni, nj);
  
}

// =========================================

// Receives a pair of coordinates and applies the grid topology to them
// Boundary conditions could be imposed here
void DaisyWorld::apply_topology(int& i, int& j) {
  
  // Toroidal topology
  if (i < 0) i = N - 1;
  else if (i >= N) i = 0;
  if (j < 0) j = N - 1;
  else if (j >= N) j = 0;

}

// =========================================

// Evaluates the noise term once
double DaisyWorld::noise_term() {
  if (noise_sigma == 0) return 0;
  else return (2.0*random_real()-1.0) * noise_sigma;
}

// =========================================

// Write a report of the current iteration to the given file
void DaisyWorld::write_report(FILE* outfile) {
  fprintf(outfile, "%d %f %f %f\n", it, global_temperature-T0, global_albedo, global_vegetation);
}

// =========================================

// Writes the current data matrices of the simulation to the given files
void DaisyWorld::writeGrids(FILE* albedo_file, FILE* temperature_file, FILE* vegetation_file) {

  int i, j, value;

  // Albedo
  fprintf(albedo_file, "# Output = %d\n", it);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {    
      if (j != 0) fprintf(albedo_file, " ");
      fprintf(albedo_file, "%f", albedo[i][j]);
    }
    fprintf(albedo_file, "\n");
  }

  // Temperature
  fprintf(temperature_file, "# Output = %d\n", it);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {    
      if (j != 0) fprintf(temperature_file, " ");
      fprintf(temperature_file, "%f", temperature[i][j]-T0);
    }
    fprintf(temperature_file, "\n");
  }

  // Vegetation
  fprintf(vegetation_file, "# Output = %d\n", it);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {    
      if (j != 0) fprintf(vegetation_file, " ");
      if (vegetation[i][j]) value = 1;
      else value = 0;
      fprintf(vegetation_file, "%d", value);
    }
    fprintf(vegetation_file, "\n");
  }


}

// =========================================
