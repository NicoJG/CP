import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

# read in the data
x,y = np.genfromtxt("exercise2.csv",delimiter=',',unpack=True,skip_header=1)
m,n = np.genfromtxt("build/exercise2_mn.csv",delimiter=',',unpack=True,skip_header=1)

# funktion for the linear regression
f = lambda x,m,n : m*x + n

# python curve fit for comparison
params,pcov = curve_fit(f,x,y)

print(f"C++ linear regression: \t m = {m}; n = {n}")
print(f"Python curve_fit: \t m = {params[0]:.6f}; n = {params[1]:.5f}")

plt.xlabel("$x$")
plt.ylabel("$y$")

plt.plot(x,f(x,m,n),'b-', label="linear regression")
plt.plot(x,y,'ko', label="data points")

plt.tight_layout()
plt.legend()
plt.savefig("build/plot_exercise2.pdf")