import sys
import os
import importlib

def load_file(path):
	if not os.path.isfile(path):
		print('Parameters file not found.')
		exit()
	module_name = os.path.basename(path).replace('-','_')
	spec = importlib.util.spec_from_loader(module_name,importlib.machinery.SourceFileLoader(module_name, path))
	module = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(module)
	sys.modules[module_name] = module
	check_file(module)
	variables = [module.r0, module.phi0, module.tmax, module.steps, module.M, module.m, 
	module.E, module.j, module.s, module.User_tensor, module.User_metric]
	if module.User_tensor:
		user_tensor = [module.g00,module.g11,module.g33,module.g_up_00,module.g_up_11,
		module.g_up_33,module.R_3001,module.R_3013,module.R_3003,module.R_3113,
		module.R_3101,module.R_1001]
		for i in user_tensor:
			variables.append(i)
	if module.User_metric:
			variables.append(module.Create_User_Metric_Tensor)

	return variables

def check_file(data):
	variables={'r0','phi0','tmax','steps','M','m','E','j','s','User_tensor','User_metric'}
	if not variables.issubset(dir(data)):
		for i in variables:
			if not {i}.issubset(dir(data)):
				print(i,' variable is missing.')
		file_template()
		exit()

	if data.User_tensor:
		user_tensor={'R_3001','R_3013','R_3003','R_3113','R_3101','R_1001','g00','g11','g33','g_up_00','g_up_11','g_up_33'}
		if not user_tensor.issubset(dir(data)):
			for i in user_tensor:
				if not {i}.issubset(dir(data)):
					print(i,' variable is missing.')
			file_template()
			exit()

	if data.User_metric and not data.User_tensor:
		user_metric={'Create_User_Metric_Tensor'}
		if not user_metric.issubset(dir(data)):
			for i in user_metric:
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
