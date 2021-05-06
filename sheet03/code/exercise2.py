import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({"font.size": 16})


eps, y = np.genfromtxt("build/ex02_energy.log", unpack=True)
name = np.arange(-20, 22, 2)
df = pd.read_csv("build/ex02_groundstate.log", header=None, sep="\s+", names=name)
df = df.pow(2)
x = np.arange(1,51)

plt.figure(figsize=(10,8))
plt.plot(eps, y, label=r"ground state energy")
plt.plot(0, -2, "r+", label="$E_{0, \mathrm{max}}$")
plt.xlabel("$\epsilon$")
plt.ylabel("$E_0(\epsilon)$")

plt.grid(ls="dotted")
plt.legend()

plt.savefig("build/ex01_gse.pdf")

plt.figure(figsize=(10,8))
plt.plot(x, df[-20], "-", label="$\epsilon = -20$")
plt.plot(x, df[0], "-", label="$\epsilon = 0$")
plt.plot(x, df[20], "-", label="$\epsilon = 20$")

plt.xlabel("$i$")
plt.ylabel(r"$| \langle \Psi_0 | i  \rangle |^2$")

plt.grid()
plt.legend()

plt.savefig("build/ex01_pd.pdf")
