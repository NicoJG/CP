import numpy as np
import matplotlib.pyplot as plt 

data = np.genfromtxt("build/exercise3_timing_data.csv",delimiter=',')

N = data[:,0]
t = data[:,1:]

for i in range(3):
    plt.plot(N,t[:,i],'-',label=f"a.{i})")

plt.xlabel("$N$")
plt.ylabel("$t \:/\: s$")

plt.legend()
plt.tight_layout()
plt.savefig("build/plot_exercise3.pdf")