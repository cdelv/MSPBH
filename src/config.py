import sys
import os
import importlib
import sympy

coord = sympy.symbols('t r theta phi')
variables={'r0','phi0','tmax','steps','M','m','E','j','s','User_tensor','User_metric','Plot_trayectory','Plot_energy'}
user_tensorv={'R_3001','R_3013','R_3003','R_3113','R_3101','R_1001','g00','g11','g33','g_up_00','g_up_11','g_up_33'}
user_metricv={'Create_User_Metric_Tensor'}

def load_file(path):
	if not os.path.isfile(path):
		print('Parameters file not found.')
		exit()
	module_name = os.path.basename(path).replace('-','_')
	spec = importlib.util.spec_from_loader(module_name,
		importlib.machinery.SourceFileLoader(module_name, path))
	module = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(module)
	sys.modules[module_name] = module
	check_file(module)
	var=[module.r0, module.phi0, module.tmax, module.steps, module.M, module.m, 
	module.E, module.j, module.s, module.User_tensor, module.User_metric,
	module.Plot_trayectory, module.Plot_energy]
	if module.User_tensor:
		user_tensor=[module.g00,module.g11,module.g33,module.g_up_00,module.g_up_11,
		module.g_up_33,module.R_3001,module.R_3013,module.R_3003,module.R_3113,
		module.R_3101,module.R_1001]
		for i in user_tensor:
			var.append(i)
	if module.User_metric:
			var.append(module.Create_User_Metric_Tensor)
	log(module)

	return var

def log(data):
	options = {False:'False', True:'True'}
	with open('log.out', 'w') as f:
		f.write('Program parameters were:\n')
		f.write('r0 = '+str(data.r0)+'\n')
		f.write('phi0 = '+str(data.phi0)+'\n')
		f.write('tmax = '+str(data.tmax)+'\n')
		f.write('steps = '+str(data.steps)+'\n')
		f.write('M = '+str(data.M)+'\n')
		f.write('m = '+str(data.m)+'\n')
		f.write('E = '+str(data.E)+'\n')
		f.write('j = '+str(data.j)+'\n')
		f.write('s = '+str(data.s)+'\n')
		f.write('Plot_trayectory = '+options[data.Plot_trayectory]+'\n')
		f.write('Plot_energy = '+options[data.Plot_energy]+'\n')
		f.write('User_tensor = '+options[data.User_tensor]+'\n')
		f.write('User_metric = '+options[data.User_metric]+'\n')
		if data.User_tensor:
			f.write('g00 = '+str(sympy.simplify(data.g00))+'\n')
			f.write('g11 = '+str(sympy.simplify(data.g11))+'\n')
			f.write('g33 = '+str(sympy.simplify(data.g33))+'\n')
			f.write('g_up_00 = '+str(sympy.simplify(data.g_up_00))+'\n')
			f.write('g_up_11 = '+str(sympy.simplify(data.g_up_11))+'\n')
			f.write('g_up_33 = '+str(sympy.simplify(data.g_up_33))+'\n')
			f.write('R_3001 = '+str(sympy.simplify(data.R_3001))+'\n')
			f.write('R_3013 = '+str(sympy.simplify(data.R_3013))+'\n')
			f.write('R_3113 = '+str(sympy.simplify(data.R_3113))+'\n')
			f.write('R_1001 = '+str(sympy.simplify(data.R_1001))+'\n')
			f.close()

	print('The program parameters were outputted to log.out.')
	

def check_file(data):
	if not variables.issubset(dir(data)):
		for i in variables:
			if not {i}.issubset(dir(data)):
				print(i,' variable is missing.')
		file_template()
		exit()

	if data.User_tensor:
		if not user_tensorv.issubset(dir(data)):
			for i in user_tensorv:
				if not {i}.issubset(dir(data)):
					print(i,' variable is missing.')
			file_template()
			exit()

	if data.User_metric and not data.User_tensor:
		if not user_metricv.issubset(dir(data)):
			for i in user_metricv:
				if not {i}.issubset(dir(data)):
					print(i,' function is missing.')
			file_template()
			exit()		

def file_template():
	print('You have to pass a file with the program parameters. Take a look at template.txt that was just created.')
	print('Edit template.txt and add the parameters you want. DO NOT CHANGE VARIABLE NAMES.')
	print('In case of an error, template.txt will be overwritten. Please use a different file for your work.')


	with open('template.txt', 'w') as f:
		f.write('''import math
import sympy
from einsteinpy.symbolic import MetricTensor, predefined
#################################
#  Initial conditions
#################################
r0 = 6
phi0 = 0

#################################
#  Propagation parameters
#################################
tmax = 50            #Maximun time
steps = 100000       #Number of steps

#################################
#  Program parameters
#################################
M = 1                #Black hole mass
m = 1                #Particle mass
E = math.sqrt(8/9)   #Energy
j = math.sqrt(12)    #Angular momentum
s = 0                #Spin
Plot_trayectory = True
Plot_energy = False

#################################
#  Metric options
#################################
User_tensor = False
User_metric = False

#################################
# User_tensor (User_tensor=True)
#################################
#These 3 components must have both indices down. (Metric tensor)
g00 = '-(1-2*'+str(M)+'/r)'
g11 = '(1-2*'+str(M)+'/r)**(-1)'
g33 = 'r**2*sin(theta)**2'
#These 3 componets must have both indices up. (Metric tensor)
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

#######################################################
# User_metric (User_tensor=Fase and User_metric=True)
#######################################################
#The coordinate order is coord = sympy.symbols('t r theta phi').
#Create the metric tensor on einsteinpy notation.
#The example below is a schwarzschild metric.

def Create_User_Metric_Tensor(coord):
	list2d = [[0 for i in range(4)] for i in range(4)]
	list2d[0][0] = -(1-2*M/coord[1])
	list2d[1][1] = 1 / (1 - 2*M / coord[1])
	list2d[2][2] = coord[1]**2
	list2d[3][3] = coord[1]**2*sympy.sin(coord[2])**2
	Metric = MetricTensor(list2d, coord)
	Metric.tensor()
	return Metric ''')
	f.close()

if __name__ == '__main__':
	main()
