#Modules
import sys
import numpy as np
import sympy
from scipy.integrate import odeint

#Source code
import src.config as config
import src.output as output
import src.propagate as propagate

#Symbolic coordinates
coord = sympy.symbols('t r theta phi')
params = sympy.symbols('t r phi')

def main():
	if len(sys.argv) != 2:
		config.file_template()
		sys.exit("Usage: python MSPBH.py template.txt")
	
	data = config.load_file(sys.argv[1])
	#Initial conditions [r, phi]
	z0=[data[0],data[1]]
	tmax = data[2]
	steps = data[3]
	#Black hole mass
	M = data[4]
	#Mass of the particle
	m = data[5] 
	#Energy of the particle
	E = data[6]
	#Orbital angular momentum
	j = data[7]
	#Spin Angular momentum
	s = data[8]

	#If user provides components of Riemann and metric tensors
	if data[9]:
		#We are always in the plain. Evaluate theta = pi/2 to speed up the code
		g00 = sympy.lambdify(params,sympy.simplify(data[13]).subs(coord[2], np.pi/2),'numpy')
		g11 = sympy.lambdify(params,sympy.simplify(data[14]).subs(coord[2], np.pi/2),'numpy')
		g33 = sympy.lambdify(params,sympy.simplify(data[15]).subs(coord[2], np.pi/2),'numpy')
		g_up_00 = sympy.lambdify(params,sympy.simplify(data[16]).subs(coord[2], np.pi/2),'numpy')
		g_up_11 = sympy.lambdify(params,sympy.simplify(data[17]).subs(coord[2], np.pi/2),'numpy')
		g_up_33 = sympy.lambdify(params,sympy.simplify(data[18]).subs(coord[2], np.pi/2),'numpy')
		R_3001 = sympy.lambdify(params,sympy.simplify(data[19]).subs(coord[2], np.pi/2),'numpy')
		R_3013 = sympy.lambdify(params,sympy.simplify(data[20]).subs(coord[2], np.pi/2),'numpy')
		R_3003 = sympy.lambdify(params,sympy.simplify(data[21]).subs(coord[2], np.pi/2),'numpy')
		R_3113 = sympy.lambdify(params,sympy.simplify(data[22]).subs(coord[2], np.pi/2),'numpy')
		R_3101 = sympy.lambdify(params,sympy.simplify(data[23]).subs(coord[2], np.pi/2),'numpy')
		R_1001 = sympy.lambdify(params,sympy.simplify(data[24]).subs(coord[2], np.pi/2),'numpy')

		#Calculate analitic partial derivatives w.r.t r
		Part_g00 = sympy.simplify(sympy.diff(sympy.simplify(data[13]).subs(coord[2], np.pi/2), coord[1]))
		Part_g33 = sympy.simplify(sympy.diff(sympy.simplify(data[15]).subs(coord[2], np.pi/2), coord[1]))
		gamma, Part_gamma = propagate.calc_gamma(sympy.simplify(data[13]).subs(coord[2], np.pi/2),
			sympy.simplify(data[14]).subs(coord[2], np.pi/2),
			sympy.simplify(data[15]).subs(coord[2], np.pi/2))

	#If user asks the program to calculate the Riemann tensor. 
	else:
		if data[10]: #User custom metric tensor
			Metric_Tensor=data[13](coord).change_config('ll')
		else: #Predefined metric tensor
			Metric_Tensor=propagate.Create_Metric_Tensor(M)

		up_Metric_Tensor = Metric_Tensor.change_config('uu') #Metric tensor with 2 up indices.
	
		#Create the Riemann curvature tensor from the metric
		Riemann_Tensor = propagate.Create_Riemann_Tensor(Metric_Tensor) #Make Riemann type (0,4)

		with open('log.out', 'a') as f:
			f.write(str(Metric_Tensor)+'\n')
			f.write(str(up_Metric_Tensor)+'\n')
			f.write(str(Riemann_Tensor)+'\n')
			f.close()


		#We are always in the plain. Evaluate theta = pi/2 to speed up the code
		g00 = sympy.lambdify(params, Metric_Tensor[0,0].subs(coord[2], np.pi/2),'numpy')
		g11 = sympy.lambdify(params, Metric_Tensor[1,1].subs(coord[2], np.pi/2),'numpy')
		g33 = sympy.lambdify(params, Metric_Tensor[3,3].subs(coord[2], np.pi/2),'numpy')
		g_up_00 = sympy.lambdify(params, up_Metric_Tensor[0,0].subs(coord[2], np.pi/2),'numpy')
		g_up_11 = sympy.lambdify(params, up_Metric_Tensor[1,1].subs(coord[2], np.pi/2),'numpy')
		g_up_33 = sympy.lambdify(params, up_Metric_Tensor[3,3].subs(coord[2], np.pi/2),'numpy')
		R_3001 = sympy.lambdify(params, Riemann_Tensor[3,0,0,1].subs(coord[2], np.pi/2),'numpy')
		R_3013 = sympy.lambdify(params, Riemann_Tensor[3,0,1,3].subs(coord[2], np.pi/2),'numpy')
		R_3003 = sympy.lambdify(params, Riemann_Tensor[3,0,0,3].subs(coord[2], np.pi/2),'numpy')
		R_3113 = sympy.lambdify(params, Riemann_Tensor[3,1,1,3].subs(coord[2], np.pi/2),'numpy')
		R_3101 = sympy.lambdify(params, Riemann_Tensor[3,1,0,1].subs(coord[2], np.pi/2),'numpy')
		R_1001 = sympy.lambdify(params, Riemann_Tensor[1,0,0,1].subs(coord[2], np.pi/2),'numpy')
		
		#Calculate analitic partial derivatives w.r.t r
		Part_g00 = sympy.simplify(sympy.diff(Metric_Tensor[0,0].subs(coord[2], np.pi/2), coord[1]))
		Part_g33 = sympy.simplify(sympy.diff(Metric_Tensor[3,3].subs(coord[2], np.pi/2), coord[1]))
		gamma, Part_gamma = propagate.calc_gamma(Metric_Tensor[0,0].subs(coord[2], np.pi/2),
			Metric_Tensor[1,1].subs(coord[2], np.pi/2),
			Metric_Tensor[3,3].subs(coord[2], np.pi/2))

	#Make derivatives and gamma fast evaluatable functions
	Part_g00 = sympy.lambdify(params,Part_g00,'numpy')
	Part_g33 = sympy.lambdify(params,Part_g33,'numpy')
	gamma = sympy.lambdify(params,gamma,'numpy')
	Part_gamma = sympy.lambdify(params,Part_gamma,'numpy')

	#Total Angular momentum
	J = j + s

	#Time span
	t=np.linspace(0,tmax,steps) #start, finish, n of points

	#Solve ODE
	z = odeint(propagate.dif_equation, z0, t, args=(m, E, s, J, gamma, Part_gamma, Part_g00, Part_g33, g00, g11, g33,
		g_up_00, g_up_11, g_up_33, R_3001, R_3013, R_3003, R_3113, R_3101, R_1001), rtol=1E-8, atol=1E-8)

	e, j = propagate.Calculate_Energy(z,t,m,E,s,J,gamma,Part_g00,Part_g33,g_up_00,g_up_11,g_up_33)
	
	output.save_data(t,z,e,j)

	#Plot trayectory and energy
	if data[11]:
		output.plot_trayectory(z,t)

	if data[12]:
		output.plot_energy(t,e,j)

if __name__ == '__main__':
	main()
