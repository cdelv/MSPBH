# MSPBH (Motion of Spinning Particles around Black Holes)

This code calculates the trajectory of a spinning particle around a black hole with an arbitrary metric tensor. The user can use the configuration file to adapt the code for a specific situation.

# Dependencies

The code uses symbolic algebra to calculate necessary quantities, such as the Riemann tensor and partial derivatives. Also uses `Odeint` to integrate equations of motion. To run this code, you need `Sympy`, `Scipy`, `Einsteinpy`, `Numpy`, and `Matplotlib`.  Use 
```
pip install -r requirements.txt
```
to install the dependencies. Make sure you are using `Python3`.

# Usage

If you execute the code like this
```
python3 MSPBH.py
```
The execution will stop, and you'll see a message with instructions. At the same time, a file called `template.txt` will be created. If the file already exists, it will be overwritten. On `template.txt`, you can find all the variables that the program needs. Make sure to use a copy of that file, otherwise, you can lose your configuration. If the variables `User_tensor` and `User_metric` are set to `False`, the program will use a premade Schwarzschild metric from the library `Einsteinpy`. If `User_tensor` is `True`, the program will use the components of the tensors given by the user on the configuration file. This option has precedence over `User_metric`. If `User_metric` is `True` and `User_tensor` is `False`, the program will use the function `Create_User_Metric_Tensor` to create a metric tensor with the library `Einsteinpy`. Then, it will compute the Riemann tensor and perform the calculations. The extension of the file doesn't matter.

When you finish editing your configuration file, type 
```
python3 MSPBH.py file_name.txt
```
The program will check for any missing variables and notify them. If this happens, `template.txt` will be overwritten. If `User_metric` is `False`, you don't have to keep the Riemann and metric tensor components of the template. They are only necessary when that option is activated. The same happens with the function `Create_User_Metric_Tensor`.

If you manage to execute the program successfully, you will see 2 output files, log.out and data.csv. The first one has information about all the variables read from the configuration file. The second one is a CSV file where you will find the time, <img src="https://render.githubusercontent.com/render/math?math=r"> coordinate, <img src="https://render.githubusercontent.com/render/math?math=\phi"> coordinate, energy, and angular momentum of the particle. This calculation is performed assuming that the particle is on the equatorial plane of the black hole (<img src="https://render.githubusercontent.com/render/math?math=\theta = \pi/2">). 

# Reference

For a complete explanation and derivation of the equations of motion, check the book "A Guide to Black Holes" from Nova Science Publishers (which hasn't been published yet).

When the book is published, I'll leave the reference here.
