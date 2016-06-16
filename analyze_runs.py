# -*- coding: utf-8 -*-
import tools
import matplotlib.pyplot as plt

# =========================

# Luminosity range
min_lumin = 0.6;
max_lumin = 1.4;
lumin_step = 0.01;
runs_per_lumin = 10;

# Simulation parameters
N = 100;
num_its = 1000;

datadir = "/home/meithan/Desktop/Daisyworld/"

# =========================  

lumins = []
avg_temps = []
teq_temps = []
  
lumin = min_lumin
while lumin <= max_lumin + lumin_step:

  lumins.append(lumin)
  teq_temps.append((lumin*917.0*(1-0.5)/ 5.670367e-8)**0.25-273.15)

  plt.clf()
  power_runs_avg = []
 
  foo = [] 
  for run in range(1,runs_per_lumin+1):
  
    fname = "L%.2f_r%02d.dat" % (lumin, run)
    print "Reading %s ..." % fname
    its, temp, albedo, veget = tools.load_data(datadir + fname)
    avg_temp = sum(temp)/float(len(temp))

    foo.append(avg_temp)
    
    # compute and plot power spectrum
    freqs, power = tools.power_spectrum(temp, sampling=1.0)
    freqs = freqs[1:]
    power = power[1:]
    if len(power_runs_avg) == 0:
      power_runs_avg = power
    else:
      for i in range(len(power_runs_avg)):
        power_runs_avg[i] += power[i]
    plt.loglog(freqs, power, color="blue", alpha=0.25)

  # Mean time-averaged temperature for the runs
  avg_temps.append(sum(foo)/float(len(foo)))

  # Compute power-law fit of average spectrum
  for i in range(len(power_runs_avg)):
    power_runs_avg[i] /= float(runs_per_lumin)
  y0, beta = tools.powerlaw_fit(freqs, power_runs_avg)

  # Finish plot for this luminosity
  plt.loglog(freqs, power_runs_avg, lw=1.5, color="black")
  plt.loglog(freqs, map(lambda f: y0*f**beta, freqs), "r--")
  plt.text(0.02, 0.05, u"β = %f" % beta, transform=plt.gca().transAxes, color="red")
  plt.grid()
  plt.title("L = %f, %i runs, 10000 iterations per run" % (lumin, runs_per_lumin))
  fname2 = fname.split("_")[0] + "_spectrum.png"
  plt.savefig(fname2)
  print "Wrote %s" % fname2
#  plt.show()
    
  lumin += lumin_step

# Make Lovelock plot

plt.clf()
plt.plot(lumins, teq_temps, "k--", lw=1.5, label="Barren world")
plt.plot(lumins, avg_temps, "b-", lw=1.5, label="Daisyworld")
plt.plot((min_lumin,max_lumin),(22.5,22.5),"b--")
ax = plt.gca()
ax.set_xlim((min_lumin,max_lumin))
plt.xlabel("Solar luminosity")
plt.ylabel(u"Planetary Temperature [°C]")
plt.legend(loc="upper left")
plt.grid()

plt.show()

