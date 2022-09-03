import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

def import_data(file_name, v):
    for l in open(file_name).readlines():
        data = l[:-1].split(' ')
        v += [float(data[0])]

# plot setting
plt.rcParams['font.family'] = 'arial'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['axes.linewidth'] = 1.2
#plt.rcParams['axes.grid'] = False
#plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.linewidth'] = 0.3
plt.rcParams["legend.markerscale"] = 2
plt.rcParams["legend.fancybox"] = False
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.edgecolor"] = 'black'

# for plot
f = plt.figure(figsize=(12.0, 8.0))
ax = f.add_subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

#ax1.axis([-1.2, 1.2, -1.2, 1.2])
#ax2.axis([0, 100, -2, 2]) 

ims = []

for i in range(0,100):
    plt.cla()
    v,v2,v3 = [],[],[]
    file_name = "fluid/" + str(i) + ".dat"
    import_data(file_name, v)
    file_name = "solid/" + str(i) + ".dat"
    import_data(file_name, v2) 
    file_name = "test/" + str(i) + ".dat"
    import_data(file_name, v3)
    plt.plot(v, color="red")   
    plt.plot(v2, color="blue")
    plt.plot(v3, color="green")   
    img_name = "img/image"+str(i)+".png"
    plt.savefig(img_name)