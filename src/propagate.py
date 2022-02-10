import numpy as np
from einsteinpy.symbolic import MetricTensor, RiemannCurvatureTensor, predefined
import sympy
coord = sympy.symbols('t r theta phi')

def Calculate_Energy(z,t,m,E,s,J,gamma,Part_g00,Part_g33,g_up_00,g_up_11,g_up_33):
	r=z[:,0]
	phi=z[:,1]
	### Evaluate gamma
	gamma=gamma(t,r,phi)
	### Evaluate metric tensor partial derivatives
	Part_g00=Part_g00(t,r,phi)
	Part_g33=Part_g33(t,r,phi)
	### Evaluate up metric tensor components
	g_up_00=g_up_00(t,r,phi)
	g_up_11=g_up_11(t,r,phi)
	g_up_33=g_up_33(t,r,phi)
	### Calculate momentum
	P0=-(2*m*gamma*(2*m*gamma*E+s*J*Part_g00))/(4*m**2*gamma**2+s**2*Part_g00*Part_g33)
	P3=(2*m*gamma*(2*m*gamma*J+s*E*Part_g33))/(4*m**2*gamma**2+s**2*Part_g00*Part_g33)
	### Calculate total energy and angular momentum
	e=-P0-s*P3*Part_g00/(2*m*gamma)
	j=P3-s*P0*Part_g33/(2*m*gamma)
	print('Energy difference is ',abs(e[0]-e[-1]))
	print('Angular momentum difference is ',abs(j[0]-j[1]))
	return e, j

def Create_Metric_Tensor(M=10):
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

	if P1_2 < 0.0000000001:
		return P0, 0.0, P3

	return P0, np.sqrt(P1_2), P3
	
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
	fac2 = s/(m*gamma**2)

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
	fac2 = s/(m*gamma**2)

	derivative = (P[2]-fac2*P[1]*drdt*Part_gamma - fac1*prod1)/(P[0]+fac1*prod2)

	return derivative

def dif_equation(z,t,m, E, s, J, gamma, Part_gamma, Part_g00, Part_g33, g00, g11, g33, g_up_00, g_up_11, g_up_33, R_3001, R_3013, R_3003, R_3113, R_3101, R_1001):
	r=z[0]
	phi=z[1]
	### Evaluate gamma and gamma derivative
	gamma = gamma(t,r,phi)
	Part_gamma = Part_gamma(t,r,phi)
	### Evaluate metric tensor partial derivatives
	Part_g00 = Part_g00(t,r,phi)
	Part_g33 = Part_g33(t,r,phi)
	### Evaluate metric tensor components
	g00 = g00(t,r,phi)
	g11 = g11(t,r,phi)
	g33 = g33(t,r,phi)
	### Evaluate up metric tensor components
	g_up_00 = g_up_00(t,r,phi)
	g_up_11 = g_up_11(t,r,phi)
	g_up_33 = g_up_33(t,r,phi)
	### Evaluate Riemann tensor Slices
	R_3001 = R_3001(t,r,phi)
	R_3013 = R_3013(t,r,phi)
	R_3003 = R_3003(t,r,phi)
	R_3113 = R_3113(t,r,phi)
	R_3101 = R_3101(t,r,phi)
	R_1001 = R_1001(t,r,phi)
	### Calculate momentum
	P = calc_momentum(g_up_00, g_up_11, g_up_33, gamma, Part_g00, Part_g33, m, E, s, J)
	### Calculate non_zero components of the spin tensor
	S = calc_spin_tensor(gamma, P, m, s)
	### Calculate derivatives 
	drdt = calc_drdt(R_3001, R_3013, R_3003, R_3113, R_3101, S, P, s, m, gamma, Part_gamma)            # r
	dphidt = calc_dphidt(R_3001, R_3013, R_3113, R_3101, R_1001, S, P, s, m, gamma, Part_gamma, drdt)  #phi

	return [drdt, dphidt]