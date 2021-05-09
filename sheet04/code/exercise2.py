import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size":16})

#load data
x, phi_1 = np.genfromtxt("build/ex02_a.log", unpack=True)
_, phi_2 = np.genfromtxt("build/ex02_b.log", unpack=True)


#plot(a)
plt.figure(figsize=(10, 8))

plt.plot(x[0:12], phi_1[0:12], "r+", label="$|x|<a$")
plt.plot(x[12:], phi_1[12:], "b+",  label="$|x|>a$")
plt.plot(x[6:], 8/x[6:], "--", label=r"$\phi_{\mathrm{mon}}$")

plt.xlabel(r"$\frac{x}{a}$")
plt.ylabel(r"$\phi_1(\frac{x}{a})$")
plt.ylim([0,10])

plt.grid(ls="dotted")
plt.legend()

plt.savefig("build/ex02_a.pdf")

#plot (b)
plt.figure(figsize=(10, 8))

plt.plot(x[0:12], phi_2[0:12], "r+", label="$|x|<a$")
plt.plot(x[12:], phi_2[12:], "b+",  label="$|x|>a$")
plt.plot(x[6:], 8/3/x[6:]**2, "--", label=r"$\phi_{\mathrm{dip}}$")

plt.xlabel(r"$\frac{x}{a}$")
plt.ylabel(r"$\phi_2(\frac{x}{a})$")
plt.ylim([0,2])

plt.grid(ls="dotted")
plt.legend()

plt.savefig("build/ex02_b.pdf")
