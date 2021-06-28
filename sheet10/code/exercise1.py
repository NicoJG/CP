import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.animation as animation

plt.rcParams.update({"font.size":16, "lines.markersize":10})

#plotting starting position
R_init = pd.read_csv("build/anim.b)set.dat", delim_whitespace=True, header=None)

plt.figure(figsize=(12,10))
plt.plot(R_init.to_numpy()[0], R_init.to_numpy()[1], "o")

plt.xlim([0,8])
plt.ylim([0,8])

plt.xlabel("$x$")
plt.ylabel("$y$")


plt.savefig("build/b)R_init.pdf") 

 
def plot(pathset, pathg, pathr, savename, h):
    t, vS, Ekin, Epot, T = np.genfromtxt("build/" +pathset, unpack=True)

    g = np.genfromtxt("build/" +pathg)
    r = np.genfromtxt("build/" + pathr)

    Ekin = np.cumsum(Ekin)/t
    Epot = np.cumsum(Epot)/t
    Etot = Ekin + Epot
    
    T = Ekin/15
    
    t *= h

    fig, axs = plt.subplots(2, 2,figsize=(12,10))
        
    axs[0, 0].plot(t, Ekin, "-", label="$E_{\mathrm{kin}}$")
    axs[0, 0].plot(t, Epot, "-", label="$E_{\mathrm{pot}}$")
    axs[0, 0].plot(t, Etot, "-", label="$E_{\mathrm{tot}}$")
    axs[0, 0].set_xlabel("$t\,[\mathrm{s}]$")
    axs[0, 0].set_ylabel("$E$")
    axs[0, 0].legend(loc="best")
    axs[0, 0].grid(ls="dotted")
    
    axs[0, 1].plot(t, T, "-", label="$T$")
    axs[0, 1].set_xlabel("$t\,[\mathrm{s}]$")
    axs[0, 1].set_ylabel("$T$")
    axs[0, 1].legend()
    axs[0, 1].grid(ls="dotted")
    
    axs[1, 0].plot(t, vS, "-", label="$v_S$")
    axs[1, 0].set_xlabel("$t\,[\mathrm{s}]$")
    axs[1, 0].set_ylabel("$v$")
    axs[1, 0].legend()
    axs[1, 0].grid(ls="dotted")
    
    axs[1, 1].plot(r, g, "-", label="$g$")
    axs[1, 1].set_xlabel("$r$")
    axs[1, 1].set_ylabel("$g(r)$")
    axs[1, 1].legend()
    axs[1, 1].grid(ls="dotted")

    plt.savefig("build/" + savename) 



plot("b)set.dat", "b)g.dat","b)r.dat", "b)set.pdf", 0.01)

plot("c)set0.01.dat", "c)g0.01.dat","c)r0.01.dat", "c)set0.01.pdf", 0.01)

plot("c)set1.dat", "c)g1.dat","c)r1.dat", "c)set1.pdf", 0.01)

plot("c)set100.dat", "c)g100.dat","c)r100.dat", "c)set100.pdf", 0.001)

plot("d)set.dat", "d)g.dat","d)r.dat", "d)set.pdf", 0.01)
