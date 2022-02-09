import numpy as np
from scipy.integrate import odeint, solve_ivp
import math
import matplotlib.pyplot as plt
import sys
import os
import importlib

import sympy
from sympy import sympify
from einsteinpy.symbolic import MetricTensor, RiemannCurvatureTensor, predefined

coord = sympy.symbols('t r theta phi')

def calc_momentum(g_up_00):
	if g_up_00 < 0.00000001:
		g_up_00 = 0
	return  math.sqrt(g_up_00) #revisar el valor del p1_2, no debe ser imaginario y tolerancia


g_up_00 = coord[1]**3

print(calc_momentum(g_up_00.subs(coord[1], 0.000000000001)))