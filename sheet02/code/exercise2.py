import numpy as np
import matplotlib.pyplot as plt 

data = np.genfromtxt("build/exercise3_timing_data.csv",delimiter=',')

N = data[:,0]
t = data[:,1:]

labels = ["a.1)","a.2) PartialPivLU","a.3) PartialPivLU","a.2) FullPivLU","a.3) FullPivLU"]

for i in range(t.shape[1]):
    plt.plot(N,t[:,i],'-',label=labels[i])

plt.xlabel("$N$")
plt.ylabel("$t \:/\: s$")

plt.legend()
plt.tight_layout()
plt.savefig("build/plot_exercise2.pdf")