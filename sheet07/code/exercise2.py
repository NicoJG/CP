import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size":16, "lines.marker":2})

#--- Ex 2.a1 2d
data_1 = np.genfromtxt("build/ex02_a_1.log")

labels = ["$x$", "$y$", "$z$", "$v_x$", "$v_y$", "$v_z$"]
rows, cols = 3, 1
fig, axs = plt.subplots(rows, cols, figsize=(12,8))
for r in range(rows):
    axs[r].plot(data_1[:, 0], data_1[:, r+1], "-")
    axs[r].plot(data_1[:, 0], data_1[:, r+4], "-")
    axs[r].legend([labels[r], labels[r+3]])   
plt.savefig("build/ex02_a_1_2d.pdf")

#--- Ex 2.a1 3d
fig = plt.figure(figsize=(12, 10))

ax = fig.add_subplot(projection="3d")
ax.view_init(25, 60)
ax.plot3D(data_1[:, 1], data_1[:, 2], data_1[:, 3])
ax.scatter(0, 0, 0, color="r")
ax.legend([r"$\vec{r}$", "center"])

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_zlabel(r"$z$")

plt.savefig("build/ex02_a_1_3d.pdf")


#--- Ex 2.a2
data_2 = np.genfromtxt("build/ex02_a_2.log")

fig = plt.figure(figsize=(12, 10))

ax = fig.add_subplot(projection="3d")
ax.view_init(25, 60)
ax.plot3D(data_2[:, 1], data_2[:, 2], data_2[:, 3])

ax.scatter(0, 0, 0, color="r")
ax.legend([r"$\vec{r}$", "center"])

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_zlabel(r"$z$")

plt.savefig("build/ex02_a_2.pdf")

#--- Ex 2.c
data_3 = np.genfromtxt("build/ex02_c.log")
t = data_3[:, 0]
E_kin = 0.5 * (data_3[:, 1]**2 + data_3[:, 2]**2 + data_3[:, 3]**2)
E_pot = 0.5 * (data_3[:, 4]**2 + data_3[:, 5]**2 + data_3[:, 6]**2)
E = E_kin + E_pot

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot()

ax.plot(t, E, "-")
ax.plot(t, E_kin, "-")
ax.plot(t, E_pot, "-")

ax.set_xlabel("$t$")
ax.set_ylabel("$E$")

ax.legend(["$E_{\mathcal{tot}}$", "$E_{\mathcal{kin}}$", "$E_{\mathcal{pot}}$" ])

plt.savefig("build/ex02_c.pdf")