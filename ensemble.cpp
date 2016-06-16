// Runs an ensemble of Daisyworld simulation for a range of luminosity
// values, making multiple runs per luminosity
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "DaisyWorld.h"
using namespace std;

// ==================================

// CONFIGURATION

// Fixed simulation parameters
const int N = 100;

// Data output
const char datadir[] = "./";
const bool output_grids = false;
const int output_interval = 100;

// Global vars 
const double T0 = 273.15;

// ==================================

// Command-line arguments
double min_lumin;
double max_lumin;
double lumin_step;
int runs_per_lumin;
int num_its;
double mutation_rate;
double noise_sigma;

// Global variables
int i, run, it;
double lumin;
time_t rawtime;
char run_fname[256], albedo_fname[256], temperature_fname[256], vegetation_fname[256], buffer[256];
FILE* albedo_file;
FILE* temperature_file;
FILE* vegetation_file;
FILE* run_file;

// ==================================

// Open new output files for a run
void open_files(int run) {
  time_t rawtime;
  time(&rawtime);
  if (output_grids) {

  }
  strcpy(run_fname, datadir);
  sprintf(buffer, "L%.2f_r%02d.dat", lumin, run);
  strcat(run_fname, buffer);
  run_file = fopen(run_fname, "w");
  fprintf(run_file, "# %s", ctime(&rawtime));
  fprintf(run_file, "# L = %.2f, run %02d\n", lumin, run);
  fprintf(run_file, "# it -- temperature -- albedo -- vegetation\n");
  
}

// Close open files
void close_files() {
  time_t rawtime;
  time(&rawtime);
  if (output_grids) {

  }
  fprintf(run_file, "# Finished %s", ctime(&rawtime));
  fclose(run_file);
}

// ==================================

int main(int argc, char** argv) {

  // Get params from command line arguments
  if (argc != 8) {
    printf("Incorrect number of arguments: %d!\n", argc-1);
    printf("Must provide: min_lumin max_lumin lumin_step runs_per_lumin num_its mutation_rate noise_sigma\n");
    exit(1);
  } else {
    min_lumin = strtod(argv[1], NULL);
    max_lumin = strtod(argv[2], NULL);
    lumin_step = strtod(argv[3], NULL);
    runs_per_lumin = strtod(argv[4], NULL);
    num_its = strtod(argv[5], NULL);
    mutation_rate = strtod(argv[6], NULL);
    noise_sigma = strtod(argv[7], NULL);
  }
  printf("Provided parameters:\n");
  printf(" Min lumin = %.3f\n", min_lumin);
  printf(" Max lumin = %.3f\n", max_lumin);
  printf(" Lumin step = %.3f\n", lumin_step);
  printf(" Runs per lumin = %d\n", runs_per_lumin);
  printf(" Iterations per run = %d\n", num_its);  
  printf(" Mutation rate = %.3f\n", mutation_rate);
  printf(" Noise sigma = %.3f\n", noise_sigma);
  printf("Data directory: %s\n", datadir);

  // Instantiate Daisyworld model
  DaisyWorld sim(N, 0.0, mutation_rate, noise_sigma);

  // Luminosity loop
  lumin = min_lumin;
  while (lumin < max_lumin*1.01) {

    sim.solar_lumin = lumin;
    printf("Executing %d runs for L = %f\n", runs_per_lumin, lumin);
    fflush(stdout);

    // Runs loop
    for (run = 1; run <= runs_per_lumin; run++) {

      printf("  Run %d ... ", run);
      fflush(stdout);
      sim.reset();
      open_files(run);
      sim.write_report(run_file);

      // Iterate simulation for this luminosity for num_its iterations
      for (it = 1; it <= num_its; it++) {
        sim.time_update();
        sim.write_report(run_file);
      }

      close_files();
      time(&rawtime);
      printf("complete - %s", ctime(&rawtime));

    }

    lumin += lumin_step;

  }

  printf("All runs complete.\n");
  
}
