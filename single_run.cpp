// Executes a single run of the Daisyworld simulation, with full data output
#include <stdio.h>
#include "DaisyWorld.h"

// ==================================

// CONFIGURATION

// Fixed simulation parameters
const int N = 100;

// Output data grids?
const bool output_grids = false;
const int output_interval = 100;

// ==================================

// Command-line arguments
double luminosity;
int num_its;
double mutation_rate;
double noise_sigma;

const double T0 = 273.15;
int i, run, it;
time_t rawtime;
char buffer[256];
FILE* albedo_file;
FILE* temperature_file;
FILE* vegetation_file;
FILE* run_file;

// ==================================

// Open new output files for a run
void open_files() {

  time_t rawtime;
  time(&rawtime);

  if (output_grids) {

    // Albedo
    sprintf(buffer, "L%.2f_albedo.dat", luminosity);
    albedo_file = fopen(buffer, "w");
    fprintf(albedo_file, "# %s", ctime(&rawtime));
    fprintf(albedo_file, "# L = %.2f\n", luminosity);
    fprintf(albedo_file, "# N = %d\n", N);
    fprintf(albedo_file, "# R = %f\n", mutation_rate);
    fprintf(albedo_file, "# noise = %f\n", noise_sigma);

    // Temperature
    sprintf(buffer, "L%.2f_temperature.dat", luminosity);
    temperature_file = fopen(buffer, "w");
    fprintf(temperature_file, "# %s", ctime(&rawtime));
    fprintf(temperature_file, "# L = %.2f\n", luminosity);
    fprintf(temperature_file, "# N = %d\n", N);
    fprintf(temperature_file, "# R = %f\n", mutation_rate);
    fprintf(temperature_file, "# noise = %f\n", noise_sigma);

    // Vegetation
    sprintf(buffer, "L%.2f_vegetation.dat", luminosity);
    vegetation_file = fopen(buffer, "w");
    fprintf(vegetation_file, "# %s", ctime(&rawtime));
    fprintf(vegetation_file, "# L = %.2f\n", luminosity);
    fprintf(vegetation_file, "# N = %d\n", N);
    fprintf(vegetation_file, "# R = %f\n", mutation_rate);
    fprintf(vegetation_file, "# noise = %f\n", noise_sigma);

  } 

  // Run file
  sprintf(buffer, "L%.2f.dat", luminosity);
  run_file = fopen(buffer, "w");
  fprintf(run_file, "# %s", ctime(&rawtime));
  fprintf(run_file, "# L = %.2f\n", luminosity);
  fprintf(run_file, "# N = %d\n", N);
  fprintf(run_file, "# R = %f\n", mutation_rate);
  fprintf(run_file, "# noise = %f\n", noise_sigma);  
  fprintf(run_file, "# it -- temperature -- albedo -- vegetation\n");
  
}

// Close open files
void close_files() {
  time_t rawtime;
  time(&rawtime);
  if (output_grids) {
    fprintf(albedo_file, "# Finished %s", ctime(&rawtime));
    fclose(albedo_file);
    fprintf(temperature_file, "# Finished %s", ctime(&rawtime));
    fclose(temperature_file);
    fprintf(vegetation_file, "# Finished %s", ctime(&rawtime));
    fclose(vegetation_file);
  }
  fprintf(run_file, "# Finished %s", ctime(&rawtime));
  fclose(run_file);
}

// ==================================

int main(int argc, char** argv) {

  // Get params from command line arguments
  if (argc != 4+1) {
    printf("Incorrect number of arguments: %d!\n", argc-1);
    printf("Must provide: luminosity, number of iterations, mutation rate, noise sigma\n");
    exit(1);
  } else {
    luminosity = strtod(argv[1], NULL);
    num_its = strtod(argv[2], NULL);
    mutation_rate = strtod(argv[3], NULL);
    noise_sigma = strtod(argv[4], NULL);
  }
  printf("Running simulation with parameters:\n");
  printf(" Luminosity = %.3f\n", luminosity);
  printf(" Iterations = %d\n", num_its);
  printf(" Mutation rate = %f\n", mutation_rate);
  printf(" Noise sigma = %f\n", noise_sigma);

  DaisyWorld sim(N, luminosity, mutation_rate, noise_sigma);
  sim.reset(12345);

  open_files();
  printf("0 %f %f %f\n", sim.global_temperature-sim.T0, sim.global_albedo, sim.global_vegetation);
  sim.write_report(run_file);
  if (output_grids) sim.writeGrids(albedo_file, temperature_file, vegetation_file);

  for (int it = 1; it <= num_its; it++) {
    sim.time_update();
    sim.write_report(run_file);
    if (it % 10000 == 0) 
      printf("%d %f %f %f\n", it, sim.global_temperature-sim.T0, sim.global_albedo, sim.global_vegetation);
    if (output_grids && (it % output_interval == 0)) {
      sim.writeGrids(albedo_file, temperature_file, vegetation_file);
    }
  }

  close_files();

}
