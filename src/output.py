import matplotlib.pyplot as plt
import numpy as np

def plot_trayectory(z,t):
	r = z[:,0]
	phi = z[:,1]
	x = r*np.cos(phi)
	y = r*np.sin(phi)
			
	plt.plot(x,y)
	plt.title('Particle trayectory', fontsize=18, pad=15)
	plt.xlabel('x', fontsize=15)
	plt.ylabel('y', fontsize=15)
	plt.show()

	plt.plot(t,z[:,0])
	plt.title('r vs time', fontsize=18, pad=15)
	plt.xlabel('Time', fontsize=15)
	plt.ylabel('r', fontsize=15)
	plt.show()

	plt.plot(t,z[:,1])
	plt.title('$\phi$ vs time', fontsize=18, pad=15)
	plt.xlabel('Time', fontsize=15)
	plt.ylabel('$\phi$', fontsize=15)
	plt.show()

def plot_energy(t,e,j):
	plt.plot(t,e)
	plt.title('Particle Energy', fontsize=18, pad=15)
	plt.xlabel('Time', fontsize=15)
	plt.ylabel('Energy', fontsize=15)
	plt.show()
	plt.plot(t,j)
	plt.title('Particle angular momentum', fontsize=18, pad=15)
	plt.xlabel('Time', fontsize=15)
	plt.ylabel('Momentum', fontsize=15)
	plt.show()

def save_data(t,z,e,j):
	h='time,r,phi,energy,momentum'
	data = np.c_[t, z[:,0], z[:,1], e, j]
	np.savetxt('data.csv', data, delimiter=',', header=h)