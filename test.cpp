// DaisyWorld class test
#include <stdio.h>
#include "DaisyWorld.h"

int main(int argc, char** argv) {

  const int num_its = 10;
  
  DaisyWorld sim(100, 1.0, 0.0, 0.0);
  sim.reset(12345);
  sim.sigma_noise = 0.0;
  
  printf("# L = %f\n", sim.solar_lumin);
  printf("0 %f %f %f\n", sim.global_temperature-sim.T0, sim.global_albedo, sim.global_vegetation);

  for (int it = 1; it <= num_its; it++) {
    sim.time_update();
    if (true) {
      printf("%d %f %f %f\n", it, sim.global_temperature-sim.T0, sim.global_albedo, sim.global_vegetation);
    }
  }


}
