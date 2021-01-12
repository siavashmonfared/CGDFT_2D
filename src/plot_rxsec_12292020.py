import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib.animation 
plt.rcParams.update({'font.size': 22})
plt.rcParams.update({'font.family': 'sans-serif'})

L = 620.;
a0 = 0.25;
L_res = 5.;

NSamp = 99;
#Ny = 2481;
#Nx = 2481;
Ny = 1241;
Nx = 1241;

IDP = 0;
IDS = 1;

data = np.zeros(Nx*Ny)
fig = plt.figure()
ax = plt.axes()
data = np.loadtxt(fname='rho_input')
data = data.reshape(Nx,Ny)
cmap = plt.get_cmap('jet')
cmap.set_bad('gray')
data = np.ma.masked_equal(data, 999)
im = ax.imshow(data, interpolation='nearest', origin='lower', cmap=cmap, vmin=0, vmax = 1)
im.set_data(data);
plt.colorbar(im,orientation='vertical',label='Density, $\\rho$');
#plt.title("$<\\rho>$ "+str(sim_data[x,1]))
plt.axis(aspect='equal')
frame1 = plt.gca()
frame1.axes.get_xaxis().set_visible(False)
frame1.axes.get_yaxis().set_visible(False)
frame1.axis('off')
plt.rcParams["figure.figsize"] = (20,3)
plt.savefig('pic.png',dpi=400)
    #plt.savefig('plot.png');

    
