# -*- coding: utf-8 -*-
import sys
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import tools

# ===========================

fig_save = True
fig_show = False
Tmin = 18
Tmax = 28

# ===========================

if len(sys.argv) != 5:
  print "Error! Must provide the following filenames as command-line arguments:"
  print "plot_simulation.py <run.dat> <albdo.dat> <temperature.dat> <vegetation.dat>"
  exit()
else:
  run_fname = sys.argv[1]
  albedo_fname = sys.argv[2]
  temperature_fname = sys.argv[3]
  vegetation_fname = sys.argv[4]

T0 = 273.15

# Read run file header
f = open(run_fname)
f.readline()
L = float(f.readline().split("=")[1])
N = int(f.readline().split("=")[1])
R = float(f.readline().split("=")[1])
noise = float(f.readline().split("=")[1])
print "L=", L
print "N=", N
print "R=", R
print "noise=", noise

# Determine output numbers from albedo data file
albedo_file = open(albedo_fname)
outputs = []
line = albedo_file.readline()
while line != "":
  if line.startswith("#") and "Output" in line:
    outputs.append(int(line.split("=")[1]))
  line = albedo_file.readline()
albedo_file.close()

its, temp, alb, veget = tools.load_data(run_fname)
albedo_file = open(albedo_fname)
temperature_file = open(temperature_fname)
vegetation_file = open(vegetation_fname)
for i in range(5):
  albedo_file.readline()
  temperature_file.readline()
  vegetation_file.readline()

# Vegetation custom 2-color colormap
veg_cmap = matplotlib.colors.ListedColormap(["white","darkgreen"])

print "Outputs found:", outputs
for out in outputs:

  albedo = []
  temperature = []
  vegetation = []
  albedo_file.readline()
  temperature_file.readline()
  vegetation_file.readline()
  for i in range(N):
    albedo.append(map(lambda x: float(x), albedo_file.readline().split()))
    temperature.append(map(lambda x: float(x), temperature_file.readline().split()))
    vegetation.append(map(lambda x: float(x), vegetation_file.readline().split()))

  fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(10,8))
  fig.suptitle("Iteration %03i" % out, fontsize=15)
  fig.subplots_adjust(left=0.04, right=0.94, bottom=0.04, top=0.92, hspace=0.15)
  
  im1 = ax1.imshow(temperature, cmap=plt.get_cmap('jet'), interpolation="none", vmin=Tmin, vmax=Tmax)
  ax1.set_title(u"Temperature [°C]")
  divider1 = make_axes_locatable(ax1)
  cax1 = divider1.append_axes("right", size="10%", pad=0.05)
  bar1 = plt.colorbar(im1, cax=cax1)

  ax2.plot(its, temp, "r-")
  ax2.plot((out),(temp[out]),"ro")
  ax2.set_ylim((Tmin,Tmax))
  ax2.set_ylabel(u"Temperature [°C]")
  ax2a = ax2.twinx()
  ax2a.plot(its, alb, color="blue")
  ax2a.plot((out), (alb[out]),"bo")
  ax2a.plot(its, veget, color="green")
  ax2a.plot((out), (veget[out]),"go")
  ax2a.set_ylim((-0.05,1.05))
  ax2a.set_ylabel("Albedo / Vegetation")
  ax2.text(0.25, 1.01, 'Temperature', va="bottom", ha="center", transform=ax2.transAxes, color="red", fontsize=14)
  ax2.text(0.515, 1.01, 'Albedo', va="bottom", ha="center", transform=ax2.transAxes, color="blue", fontsize=14)
  ax2.text(0.75, 1.01, 'Vegetation', va="bottom", ha="center", transform=ax2.transAxes, color="green", fontsize=14)
  
  im3 = ax3.imshow(vegetation, cmap=veg_cmap, interpolation="none", vmin=0, vmax=1)
  ax3.set_title("Vegetation")
  divider3 = make_axes_locatable(ax3)
  cax3 = divider3.append_axes("right", size="10%", pad=0.05, aspect=2, anchor=(0,1))
  bar3 = fig.colorbar(im3, cax=cax3)
  bar3.set_ticks((0.25, 0.75))
  bar3.set_ticklabels(["no","yes"])

  im4 = ax4.imshow(albedo, cmap=plt.get_cmap('afmhot'), interpolation="none", vmin=0, vmax=1)
  ax4.set_title("Albedo")
  divider4 = make_axes_locatable(ax4)
  cax4 = divider4.append_axes("right", size="10%", pad=0.05)
  bar4 = plt.colorbar(im4, cax=cax4)

  if fig_save:
    fname = "plot%05i.png" % out
    plt.savefig(fname)
    print "Wrote %s" % fname
  if fig_show:
    plt.show()
  
  plt.close(fig)
