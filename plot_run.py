import matplotlib.pyplot as plt

fig_save = True
fig_show = False
Tmin = 18
Tmax = 28

f = open(sys.argv[1])
its = []
temperature = []
albedo = []
vegetation = []
for line in f:
  if line.startswith("#"): continue
  data = line.strip().split()
  its.append(int(data[0]))
  temperature.append(int(data[1]))
  albedo.append(int(data[2]))
  vegetation.append(int(data[3]))
f.close()  
  
plt.plot(its, temp, "r-")
ax.set_ylim((Tmin,Tmax))
ax.set_ylabel(u"Temperature [°C]")
ax2 = ax.twinx()
ax2.plot(its, alb, color="blue")
ax2.plot((out), (alb[out]),"bo")
ax2.plot(its, veget, color="green")
ax2.plot((out), (veget[out]),"go")
ax2.set_ylim((-0.05,1.05))
ax2.set_ylabel("Albedo / Vegetation")
ax.text(0.25, 1.01, 'Temperature', va="bottom", ha="center", transform=ax.transAxes, color="red", fontsize=14)
ax.text(0.515, 1.01, 'Albedo', va="bottom", ha="center", transform=ax.transAxes, color="blue", fontsize=14)
ax.text(0.75, 1.01, 'Vegetation', va="bottom", ha="center", transform=ax.transAxes, color="green", fontsize=14)

plt.show()
