/*========================\\
|| Daisy WorldModel Class ||
\\========================*/
#ifndef DAISYWORLD_H
#define DAISYWORLD_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

class DaisyWorld {

  public:

  /*==================================================================*/

  /* MEMBER VARIABLES */

  // Simulation static parameters
  static const double sigmaB = 5.670367e-8;    // Stefan-Boltzmann constant, W m^-2 K^-4
  static const double S0 = 917.0;        // Solar flux at L=1.0, W/m^2; Watson & Lovelock value, notice this is less than Earth
  static const double T0 = 273.15;       // 0°C in Kelvin
  static const double DT = 500.0;        // Temperature diffusion constant
  static const double CT = 2500.0;       // Temperature dynamics scaling
  static const double A0 = 0.5;          // Albedo of bare ground
  static const double gamma = 0.02;      // Daisy death probability
  static const double Tmin = 5.0+273.15;     // Minimum growth temperature, K
  static const double Tmax = 40.0+273.15;    // Maximum growth temperature, K
  static const double Tideal = 22.5+273.15;  // Ideal growth growth temperature, K

  // Constructor-defined simulation parameters
  int N;                 // Dimension of the grid (N x N)
  double mutation_rate;  // Mutation amplitude
  double noise_sigma;    // Noise amplitude  
  double solar_lumin;    // Solar luminosity factor

  // Other parameters
  int num_cells;            // Number of cells (=N^2)
  double growth_constant;   // Constant for the growth law
  double Teq;               // Thermodynamic equilibrium temperature

  // RNG stuff
  int global_seed;           // Global seed  

  // Simulation variables
  int it;
  double global_temperature;
  double global_albedo;
  double global_vegetation;

  // Indicates whether cell (i,j) is covered in vegetation (daisies) or not
  bool** vegetation;
  bool** vegetation_buffer;

  // Local cell temperature (and a buffer for temperature dynamics)
  double** temperature;
  double** temperature_buffer;

  // Local cell albedo (and a buffer for albedo dynamics)
  double** albedo;
  double** albedo_buffer;
  
  /*==================================================================*/

  /* MEMBER FUNCTIONS */
  DaisyWorld(int N, double solar_lumin, double R, double sigma_noise);
  void reset();
  void reset(int seed);
  double random_real();
  double noise_term();
  double compute_global_temperature();
  double compute_global_albedo();
  double compute_global_vegetation();
  void time_update();
  void update_vegetation();
  void update_temperature();
  double child_albedo(double prog_albedo);
  double compute_beta(double T);
  void get_neigh_coords(int i, int j, int neigh, int& ni, int& nj);
  void apply_topology(int& i, int& j);
  void write_report(FILE* outfile);
  void writeGrids(FILE* albedo_file, FILE* temperature_file, FILE* vegetation_file);
};

#endif // DAISYWORLD_H
