import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np

def plot_trayectory(z,t):
	r = z[:,0]
	phi = z[:,1]
	x = r*np.cos(phi)
	y = r*np.sin(phi)

	fig, axs = plt.subplot_mosaic([['a)', 'b)'],['a)', 'c)']],figsize=(11,5),constrained_layout=True)
	axs['a)'].set_title('Particle Trajectory', fontsize=18, pad=15,fontstyle='italic')
	circle1 = plt.Circle((0, 0), r[0]/17, color='orange')
	circle2 = plt.Circle((0, 0), r[0]/18, color='red',alpha=0.7)
	circle3 = plt.Circle((0, 0), r[0]/20, color='black')
	axs['a)'].add_patch(circle1)
	axs['a)'].add_patch(circle2)
	axs['a)'].add_patch(circle3)
	axs['a)'].set_ylabel('$y$')
	axs['a)'].set_xlabel('$x$')
	axs['a)'].plot(x, y,color='blue',linewidth=2.5, alpha=0.65)
	axs['a)'].grid(color='gray', linestyle='--',linewidth=0.7,alpha=0.5)
	axs['b)'].set_title('$r$ vs. Time', fontsize=18, pad=15,fontstyle='italic')
	axs['b)'].set_ylabel('$r$')
	axs['b)'].tick_params(labelbottom=False)
	axs['b)'].plot(t, r,color='orange',linewidth=3,alpha=0.8)   
	axs['b)'].grid(color='gray', linestyle='--',linewidth=0.7,alpha=0.5) 
	axs['c)'].set_title('$\phi$ vs. Time', fontsize=18, pad=15,fontstyle='italic')
	axs['c)'].set_xlabel('Time')
	axs['c)'].set_ylabel('$\phi$')
	axs['c)'].plot(t, phi,color='green',linewidth=3,alpha=0.8)
	axs['c)'].grid(color='gray', linestyle='--',linewidth=0.7,alpha=0.5)

	for label, ax in axs.items():
	    ax.set_title(label, fontfamily='serif', loc='left', fontsize='medium')

	fig.savefig('trajectory.eps', format='eps', dpi=600)
	plt.show()

def plot_energy(t,e,j):
	fig, axs = plt.subplot_mosaic([['a)', 'b)']],figsize=(11,5),constrained_layout=True)
	axs['a)'].set_ylabel('Energy')
	axs['a)'].set_xlabel('Time')
	axs['a)'].set_title('Total Energy vs. Time', fontsize=18, pad=15,fontstyle='italic')
	axs['a)'].plot(t, e,color='orange',linewidth=2.5, alpha=0.8)
	axs['a)'].grid(color='gray', linestyle='--',linewidth=0.7,alpha=0.5)
	axs['b)'].set_ylabel('Momentum')
	axs['b)'].set_xlabel('Time')
	axs['b)'].set_title('Total Angular Momentum vs. Time', fontsize=18, pad=15,fontstyle='italic')
	axs['b)'].plot(t, j,color='red',linewidth=2.5, alpha=0.8)
	axs['b)'].grid(color='gray', linestyle='--',linewidth=0.7,alpha=0.5)

	for label, ax in axs.items():
	    ax.set_title(label, fontfamily='serif', loc='left', fontsize='medium')

	fig.savefig('energy.eps', format='eps', dpi=600)
	plt.show()

def save_data(t,z,e,j):
	h='time,r,phi,energy,momentum'
	data = np.c_[t, z[:,0], z[:,1], e, j]
	np.savetxt('data.csv', data, delimiter=',', header=h)