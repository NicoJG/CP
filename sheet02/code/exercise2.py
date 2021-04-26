import numpy as np
import matplotlib.pyplot as plt

t1, t2, t3 = np.genfromtxt('times.dat',unpack=True)
x = np.linspace(1,1000, 999)


plt.yscale('log')
plt.xscale('log')
plt.xlabel('Size N')
plt.ylabel('time[s]')
plt.plot(x, t1, "g-",label= "t1")
plt.plot(x, t2, "m-",label="t2")
plt.plot(x, t3 , "b-", label ="t3")
plt.legend(loc='best')

plt.savefig("build/ex02_plots.pdf")

