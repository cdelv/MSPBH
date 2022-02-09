import math
import sys
import os
import importlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Symbolic algebra libraries
import sympy
from sympy import sympify
from einsteinpy.symbolic import MetricTensor, RiemannCurvatureTensor, predefined

#Symbolic coordinates
coord = sympy.symbols('t r theta phi')

# To do:
# Energy plots conditions on template
# A good .out log


def main():
	if len(sys.argv) != 2:
		file_template()
		sys.exit("Usage: python acresion.py parameters.txt")
	
	data = load_file(sys.argv[1])

	#Black hole mass
	M = data[3]

	#Initial conditions   [r, phi]
	z0=[data[0],data[1]]
	tmax = data[7]
	steps = data[8]

	#Mass of the particle
	m = data[4] 

	#Energy of the particle
	E = data[2]

	#Orbital angular momentum
	j = data[5]

	#Spin Angular momentum
	s = data[6]

	#If user gives components of Riemann and metric tensors
	if data[9]:
		#We are always in the plain. Evaluate theta = pi/2 to speed up the code
		R_3001 = sympy.simplify(data[10]).subs(coord[2], math.pi/2)
		R_3013 = sympy.simplify(data[11]).subs(coord[2], math.pi/2)
		R_3003 = sympy.simplify(data[12]).subs(coord[2], math.pi/2)
		R_3113 = sympy.simplify(data[13]).subs(coord[2], math.pi/2)
		R_3101 = sympy.simplify(data[14]).subs(coord[2], math.pi/2)
		R_1001 = sympy.simplify(data[15]).subs(coord[2], math.pi/2)
		g00 = sympy.simplify(data[16]).subs(coord[2], math.pi/2)
		g11 = sympy.simplify(data[17]).subs(coord[2], math.pi/2)
		g33 = sympy.simplify(data[18]).subs(coord[2], math.pi/2)
		g_up_00 = sympy.simplify(data[19]).subs(coord[2], math.pi/2)
		g_up_11 = sympy.simplify(data[20]).subs(coord[2], math.pi/2)
		g_up_33 = sympy.simplify(data[21]).subs(coord[2], math.pi/2)

	#If user asks the program to calculate the Riemann tensor. 
	else:
		#Create the metric tensor
		if data[22]: #User custom metric tensor
			Metric_Tensor=data[23].change_config('ll')
		else: #Predefined metric tensor
			Metric_Tensor=Create_Metric_Tensor(M)

		up_Metric_Tensor = Metric_Tensor.change_config('uu') #Metric tensor with 2 up indices.

		print(Metric_Tensor)
		print(up_Metric_Tensor)
	
		#Create the Riemann curvature tensor from the metric
		Riemann_Tensor = Create_Riemann_Tensor(Metric_Tensor) #Make Riemann type (0,4)

		#We are always in the plain. Evaluate theta = pi/2 to speed up the code
		R_3001 = Riemann_Tensor[3,0,0,1].subs(coord[2], math.pi/2)
		R_3013 = Riemann_Tensor[3,0,1,3].subs(coord[2], math.pi/2)
		R_3003 = Riemann_Tensor[3,0,0,3].subs(coord[2], math.pi/2)
		R_3113 = Riemann_Tensor[3,1,1,3].subs(coord[2], math.pi/2)
		R_3101 = Riemann_Tensor[3,1,0,1].subs(coord[2], math.pi/2)
		R_1001 = Riemann_Tensor[1,0,0,1].subs(coord[2], math.pi/2)
		g00 = Metric_Tensor[0,0].subs(coord[2], math.pi/2)
		g11 = Metric_Tensor[1,1].subs(coord[2], math.pi/2)
		g33 = Metric_Tensor[3,3].subs(coord[2], math.pi/2)
		g_up_00 = up_Metric_Tensor[0,0].subs(coord[2], math.pi/2)
		g_up_11 = up_Metric_Tensor[1,1].subs(coord[2], math.pi/2)
		g_up_33 = up_Metric_Tensor[3,3].subs(coord[2], math.pi/2)
		

	#Calculate analitic Metric partial derivatives w.r.t r
	Part_g00 = sympy.simplify(sympy.diff(g00, coord[1]))
	Part_g33 = sympy.simplify(sympy.diff(g33, coord[1]))
	gamma, Part_gamma = calc_gamma(g00,g11,g33)
	#print_variables(g00, g11, g33, g_up_00, g_up_11, g_up_33, R_3001, R_3013, R_3003, R_3113, R_3101, R_1001)


	#Total Angular momentum
	J = j + s

	#Time span
	t=np.linspace(0,tmax,steps) #start, finish, n of points

	#Solve ODE
	z = odeint(dif_equation, z0, t, args=(m, E, s, J, gamma, Part_gamma, Part_g00, Part_g33, g00, g11, g33,
		g_up_00, g_up_11, g_up_33, R_3001, R_3013, R_3003, R_3113, R_3101, R_1001), rtol=1E-8, atol=1E-8)

	#Plot trayectory and energy
	plot_trayectory(z,t) #r vs phi
	plot_energy(z,t,m,E,s, J, gamma, Part_g00, Part_g33,g_up_00, g_up_11, g_up_33)

def print_variables(g00, g11, g33, g_up_00, g_up_11, g_up_33, R_3001, R_3013, R_3003, R_3113, R_3101, R_1001):
	#print('g_11=',sympy.pprint(g11))
	sympy.print(g11,"g_11")
	#sympy.pprint('g_33='+str(g33))
	#sympy.pprint('g_up_00='+str(g_up_00))
	#sympy.pprint('g_up_11='+str(g_up_11))
	#sympy.pprint('g_up_33='+str(g_up_33))
	#sympy.pprint('R_3001=',R_3001)
	#sympy.pprint('R_3013=',R_3013)
	#sympy.pprint('R_3003=',R_3003)
	#sympy.pprint('R_3113=',R_3113)
	#sympy.pprint('R_3101=',R_3101)
	#sympy.pprint('R_1001=',R_1001)

def Create_Metric_Tensor(M=10):
	#### Define problem metric #########
	#### schwarzschild metric ##########
	'''
	list2d = [[0 for i in range(4)] for i in range(4)]
	list2d[0][0] = -(1-2*M/coord[1])
	list2d[1][1] = 1 / (1 - 2*M / coord[1])
	list2d[2][2] = coord[1]**2
	list2d[3][3] = coord[1]**2*sympy.sin(coord[2])**2
	Metric = MetricTensor(list2d, coord)
	Metric.tensor()
	'''
	########## Premade schwarzschild metric example ############# c, G, a, M
	Metric = predefined.janis_newman_winicour.JanisNewmanWinicour(1, 1, 1, M)
	return Metric  ### Metric_{ab}  a, b in (0 to 4), 2 indices down

def Create_Riemann_Tensor(Metric):
	#### Calculate Riemann Tensor from the metric tensor #########
	Riemann = RiemannCurvatureTensor.from_metric(Metric)
	Riemann.tensor()
	return Riemann.change_config('llll', Metric)

def calc_gamma(g00,g11,g33):
	#Calculates the analitic expresion of gamma and its derivative w.r.t. to r
	gamma = sympy.simplify(sympy.sqrt(-1*g00*g11*g33))
	Part_gamma = sympy.simplify(sympy.diff(gamma, coord[1]))

	return gamma, Part_gamma

def calc_momentum(g_up_00, g_up_11, g_up_33, gamma, Part_g00, Part_g33, m, E, s, J):
	P0 = -(2*m*gamma*(2*m*gamma*E+s*J*Part_g00))/(4*m**2*gamma**2+s**2*Part_g00*Part_g33)
	P3 = (2*m*gamma*(2*m*gamma*J+s*E*Part_g33))/(4*m**2*gamma**2+s**2*Part_g00*Part_g33)
	P1_2 = -(m**2+g_up_00*P0**2+g_up_33*P3**2)/g_up_11

	if P1_2 < 0.000000001:
		return P0, 0.0, P3

	return P0, math.sqrt(P1_2), P3
	
def calc_spin_tensor(gamma, P, m, s):
	num = -s/(m*gamma)
	S_13 = P[0]*num
	S_01 = P[2]*num
	S_03 = P[1]*num

	return S_13, S_01, S_03

def calc_drdt(R_3001, R_3013, R_3003, R_3113, R_3101, S, P, s, m, gamma, Part_gamma):
	### Calculate R_30ab * S^ab
	prod1 = R_3013*S[0] + R_3001*S[1] + R_3003*S[2]
	### Calculate R_31ab * S^ab
	prod2 = R_3113*S[0] + R_3101*S[1] + R_3013*S[2] 
	#                                   R_3103=R_3013
	fac1 = s/(2*m*gamma)
	fac2 = s/(2*m*gamma**2)

	derivative = (P[1] + fac1*prod1)/(P[0]-P[2]*fac2*Part_gamma-fac1*prod2)

	return derivative

def calc_dphidt(R_3001, R_3013, R_3113, R_3101, R_1001, S, P, s, m, gamma, Part_gamma, drdt):
	### Calculate R_10ab * S^ab
	prod1 = R_3101*S[0] + R_1001*S[1] + R_3001*S[2]
	#       R_1013=R_3101,              R_1003=R_3001
	### Calculate R_31ab * S^ab
	prod2 = R_3113*S[0] -  R_3101*S[1]  -  R_3013*S[2]
	#       R_1313=R_3113, R_1301=-R_3101, R_1303=-R_3013     
	fac1 = s/(2*m*gamma)
	fac2 = s/(2*m*gamma**2)

	derivative = (P[2]-fac2*P[1]*drdt*Part_gamma - fac1*prod1)/(P[0]+fac1*prod2)

	return derivative

def dif_equation(z,t,m, E, s, J, gamma, Part_gamma, Part_g00, Part_g33, g00, g11, g33, g_up_00, g_up_11, g_up_33, R_3001, R_3013, R_3003, R_3113, R_3101, R_1001):
	r=z[0]
	phi=z[1]
	### Evaluate gamma and gamma derivative
	gamma = gamma.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	Part_gamma = Part_gamma.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	### Evaluate metric tensor partial derivatives
	Part_g00 = Part_g00.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	Part_g33 = Part_g33.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	### Evaluate metric tensor components
	g00 = g00.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	g11 = g11.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	g33 = g33.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	### Evaluate up metric tensor components
	g_up_00 = g_up_00.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	g_up_11 = g_up_11.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	g_up_33 = g_up_33.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	### Evaluate Riemann tensor Slices
	R_3001 = R_3001.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	R_3013 = R_3013.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	R_3003 = R_3003.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	R_3113 = R_3113.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	R_3101 = R_3101.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	R_1001 = R_1001.subs([(coord[0], t), (coord[1], r), (coord[3], phi)])
	### Calculate momentum
	P = calc_momentum(g_up_00, g_up_11, g_up_33, gamma, Part_g00, Part_g33, m, E, s, J)
	### Calculate non_zero components of the spin tensor
	S = calc_spin_tensor(gamma, P, m, s)
	### Calculate derivatives 
	drdt = calc_drdt(R_3001, R_3013, R_3003, R_3113, R_3101, S, P, s, m, gamma, Part_gamma)        # r
	dphidt = calc_dphidt(R_3001, R_3013, R_3113, R_3101, R_1001, S, P, s, m, gamma, Part_gamma, drdt)  #phi

	return [drdt, dphidt]

def plot_trayectory(z,t):
	r = z[:,0]
	phi = z[:,1]
	x = []
	y = []
	## Convert results to cartesian coordinates for the plot
	for i in range(len(r)):
		x.append(r[i]*math.cos(phi[i]))
		y.append(r[i]*math.sin(phi[i]))
			
	plt.plot(np.asarray(x),np.asarray(y))
	plt.title('Particle trayectory', fontsize=18, pad=15)
	plt.xlabel('x', fontsize=15)
	plt.ylabel('y', fontsize=15)
	plt.show()

	plt.plot(t,np.asarray(x))
	plt.title('r coordinate vs time', fontsize=18, pad=15)
	plt.xlabel('time', fontsize=15)
	plt.ylabel('r', fontsize=15)
	plt.show()

	plt.plot(t,np.asarray(y))
	plt.title('phi coordinate vs time', fontsize=18, pad=15)
	plt.xlabel('time', fontsize=15)
	plt.ylabel('phi', fontsize=15)
	plt.show()

def plot_energy(z,t,m,E,s, J, gamma, Part_g00, Part_g33,g_up_00, g_up_11, g_up_33):
	#E=-P0-s*P3*Part_g00/(2*m*gamma)
	e=np.zeros(len(t))
	j=np.zeros(len(t))
	P0=np.zeros(len(t))
	P3=np.zeros(len(t))
	PG00=np.zeros(len(t))
	PG33=np.zeros(len(t))
	Gamma=np.zeros(len(t))
	for i in range(len(t)):
		r=z[i,0]
		phi=z[i,1]
		### Evaluate gamma
		GGamma = gamma.subs([(coord[0], t[i]), (coord[1], r), (coord[3], phi)])
		### Evaluate metric tensor partial derivatives
		PPart_g00 = Part_g00.subs([(coord[0], t[i]), (coord[1], r), (coord[3], phi)])
		PPart_g33 = Part_g33.subs([(coord[0], t[i]), (coord[1], r), (coord[3], phi)])
		### Evaluate up metric tensor components
		G_up_00 = g_up_00.subs([(coord[0], t[i]), (coord[1], r), (coord[3], phi)])
		G_up_11 = g_up_11.subs([(coord[0], t[i]), (coord[1], r), (coord[3], phi)])
		G_up_33 = g_up_33.subs([(coord[0], t[i]), (coord[1], r), (coord[3], phi)])

		### Calculate momentum
		P0[i]=-(2*m*GGamma*(2*m*GGamma*E+s*J*PPart_g00))/(4*m**2*GGamma**2+s**2*PPart_g00*PPart_g33)
		P3[i]=(2*m*GGamma*(2*m*GGamma*J+s*E*PPart_g33))/(4*m**2*GGamma**2+s**2*PPart_g00*PPart_g33)
		PG00[i]=PPart_g00
		PG33[i]=PPart_g33
		Gamma[i]=GGamma

	e=-P0-s*P3*PG00/(2*m*Gamma)
	j=P3-s*P0*PG33/(2*m*GGamma)
	print('Energy difference is ', e[0]-e[len(t)-1])
	print('Angular momentum difference is ', j[0]-j[len(t)-1])

	plt.plot(t,e)
	plt.title('Particle Energy', fontsize=18, pad=15)
	plt.xlabel('t', fontsize=15)
	plt.ylabel('Energy', fontsize=15)
	plt.show()
	plt.plot(t,j)
	plt.title('Particle angular momentum', fontsize=18, pad=15)
	plt.xlabel('t', fontsize=15)
	plt.ylabel('momentum', fontsize=15)
	plt.show()


def load_file(path):
	if not os.path.isfile(path):
		print('Parameters file not found.')
		exit()

	module_name = os.path.basename(path).replace('-', '_')
	spec = importlib.util.spec_from_loader(
	    module_name,
	    importlib.machinery.SourceFileLoader(module_name, path)
	)
	module = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(module)
	sys.modules[module_name] = module
	check_file(module)
	variables = [module.r0, module.phi0, module.E, module.M, module.m, 
	module.j, module.s, module.tmax, module.steps, module.User_tensor, module.R_3001,
	module.R_3013,module.R_3003,module.R_3113,module.R_3101,module.R_1001
	,module.g00,module.g11,module.g33,module.g_up_00,module.g_up_11, module.g_up_33, module.User_metric]
	if module.User_metric and not module.User_tensor:
		variables.append(module.Create_User_Metric_Tensor(coord))
	return variables

def check_file(data):
	variables = {'r0', 'phi0', 'E', 'M','m', 'j', 's', 'tmax', 'steps', 'User_tensor', 'User_metric'}
	if not variables.issubset(dir(data)):
		for i in variables:
			if not {i}.issubset(dir(data)):
				print(i,' variable is missing.')
		file_template()
		exit()

	variables2 = {'R_3001', 'R_3013', 'R_3003', 'R_3113', 'R_3101', 'R_1001', 'g00', 'g11', 
	'g33', 'g_up_00', 'g_up_11', 'g_up_33'}
	if data.User_tensor:
		if not variables.issubset(dir(data)):
			for i in variables2:
				if not {i}.issubset(dir(data)):
					print(i,' variable is missing.')
			file_template()
			exit()

	variables3 = {'Create_User_Metric_Tensor'}
	if data.User_metric:
		if not variables3.issubset(dir(data)):
			for i in variables3:
				if not {i}.issubset(dir(data)):
					print(i,' variable is missing.')
			file_template()
			exit()

def file_template():
	print('You have to pass a file with the program parameters. Take a look at the file template.txt that was just created.')
	print('Edit template.txt and add the parameters you want. Do not change variable names.')
	print('In case of an error, template.txt will be overwritten. Please use a different name for your working file.')


	with open('template.txt', 'w') as f:
		f.write('''import math
import sympy
from einsteinpy.symbolic import MetricTensor, predefined
#######################
# Initial conditions
#######################
r0 = 6 
phi0 = 0

#######################
# Program parameters
#######################
M = 1 #Black hole mass
m = 1   #Particle mass
E = math.sqrt(8/9)   #Energy
j = math.sqrt(12)   #Angular momentum
s = 0   #Spin

#######################
# Propagation parameters
#######################
tmax = 50        #Maximun time
steps = 100000   #Number of steps

#######################
#   Metric
#######################
#If User_tensor = True, the program uses the components passed by the user. 
#If User_tensor = False and User_metric = False, the program uses predetermine einsteinpy JanisNewmanWinicour(1, 1, 1, M) metric.
#If User_tensor = True and User_metric = True, the program uses the components passed by the user.
#If User_tensor = Fase and User_metric = True, the program uses function Create_User_Metric_Tensor(coord) to calculate the Riemann tensor from the defined metric. 
User_tensor = False
User_metric = False
#######################
#   User_tensor  (User_tensor = True)
#######################
#Enter components in Sympy compatible notation. Remember only accepted coordinates are ('t r theta phi')
#Symbolic constants are not admitted.
#Metric tensor components must be in configuration [-1,1,1,1].
#Input in this section must be strings. 

#These 3 components must have both indices down.
g00 = '-(1-2*'+str(M)+'/r)'
g11 = '(1-2*'+str(M)+'/r)**(-1)'
g33 = 'r**2*sin(theta)**2'

#These 3 componets must have both indices up.
g_up_00 = '-(1-2*'+str(M)+'/r)**(-1)'
g_up_11 = '(1-2*'+str(M)+'/r)'
g_up_33 = '(r**2*sin(theta))**(-1)'

#Riemann tensor components with all 4 indices down.
R_3001 = '0'
R_3013 = '0'
R_3003 = '-(1-2*'+str(M)+'/r)*sin(theta)**2/r*'+str(M)
R_3113 = '-('+str(M)+'*sin(theta)**2)/(2*'+str(M)+'-r)'
R_3101 = '0'
R_1001 = str(M)+'*2/r**3' 

#######################
#   User_metric  (User_tensor = Fase and User_metric = True)
#######################
#Define a function from which the program will create your metric and then calculate the Riemann tensor with it.
#Use Sympy compatible notation. remember that only accepted variables are coord = sympy.symbols('t r theta phi').
#Use array notation tu acces the variables, t=coord[0], r=coord[1], theta=coord[2], phi=coord[3].
#This is an einsteinpy metric tensor.

def Create_User_Metric_Tensor(coord):
	list2d = [[0 for i in range(4)] for i in range(4)]
	list2d[0][0] = -(1-2*M/coord[1])
	list2d[1][1] = 1 / (1 - 2*M / coord[1])
	list2d[2][2] = coord[1]**2
	list2d[3][3] = coord[1]**2*sympy.sin(coord[2])**2
	Metric = MetricTensor(list2d, coord)
	Metric.tensor()
	return Metric 


#Please report any bugs to cdelv@unal.edu.co or to the git repository https://github.com/cdelv/Astrophysics''')
	f.close()


if __name__ == '__main__':
	main()
