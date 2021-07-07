import numpy as np
import matplotlib.pyplot as plt

H, m_mc = np.genfromtxt('build/mc_mag.dat',unpack=True)
m_an = np.tanh(H)

plt.xlabel('external magnetic field H')
plt.ylabel('magnetization m')
plt.plot(H, m_mc,label= "MC", alpha= 0.8)
plt.plot(H, m_an,label="analytical",alpha= 0.5)

plt.legend(loc='best')

plt.savefig("build/ex01_plot.pdf")


