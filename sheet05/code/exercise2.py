import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size":16})

for i in range(3, 5):
    x1, y1, x2, y2 = np.genfromtxt("build/ex02_m" + str(i) + ".log", unpack=True)
    k = np.arange(0, len(x1))
    plt.figure(figsize=(16,6))
    plt.subplot(121)
    plt.title("real part")
    plt.plot(k, x1, "+", label="FFT")
    plt.plot(k, x2, "x", label="FT")

    plt.xlabel("k")
    plt.ylabel(r"$\mathcal{R}(F(K))$")

    plt.grid(ls="dotted")
    plt.legend()

    plt.subplot(122)
    plt.title("imaginary part")
    plt.plot(k, y1, "+", label="FFT")
    plt.plot(k, y2, "x", label="FT")

    plt.xlabel("k")
    plt.ylabel(r"$\mathcal{Im}(F(K))$")

    plt.grid(ls="dotted")
    plt.legend()
    plt.savefig("build/ex02_m" + str(i) + ".pdf")    
