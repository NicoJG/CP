import numpy as np
import matplotlib.pyplot as plt


plt.rcParams.update({"font.size":16})

def f(x):
    return x**2 - 2
i_1, min_1 , i_2, min_2 = np.genfromtxt("build//ex03.log", unpack=True)

x = np.linspace(-1, 1, 100)


plt.figure(figsize=(10, 8))
plt.plot(x, f(x), label="$f(x)$")
plt.plot(min_1, f(min_1), "rx", ms=10, label="Bisection, $i=33$")
plt.plot(min_2, f(min_2), "b+", ms=15, label="Newton, $i=i$")

plt.xlabel("x")
plt.ylabel("f(x)")

plt.grid(ls="dotted")
plt.legend()

plt.savefig("build/ex03.pdf")    
