import numpy as np
import matplotlib.pyplot as plt 

x_g, y_g = np.genfromtxt('build/ex02gd.dat',unpack=True)
x_cg, y_cg = np.genfromtxt('build/ex02cgd.dat',unpack=True)

X, Y = np.meshgrid(x_g, y_g)
Z = (1-X)**2 + 100*(Y - X**2)**2
plt.contourf(X,Y,Z)
plt.colorbar()
plt.xlabel("x")
plt.ylabel('y')

plt.savefig("build/ex02gd.pdf") 
plt.clf()
Q, R = np.meshgrid(x_cg, y_cg)
T = (1-Q)**2 + 100*(R - Q**2)**2
cs = plt.contourf(Q,R,T)
plt.xlabel("x")
plt.ylabel('y')
plt.colorbar()

plt.savefig("build/ex02cgd.pdf")
plt.clf()
err_g = []
for i in range(len(x_g)):
    err_g.append(np.linalg.norm([x_g[i]-1,y_g[i]-1]))

x = np.linspace(0,len(err_g),len(err_g))
plt.plot(x, err_g)

plt.xlabel("i")
plt.ylabel('L^2-Norm')
plt.tight_layout()
plt.savefig("build/ex02gd_err.pdf")    
plt.clf()


err_cg = []
for i in range(len(x_cg)):
    err_cg.append(np.linalg.norm([x_cg[i]-1,y_cg[i]-1]))

r = np.linspace(0,len(err_cg),len(err_cg))
plt.plot(r, err_cg)
 
plt.xlabel("i")
plt.ylabel('L^2-Norm')  
plt.tight_layout() 
plt.savefig("build/ex02cgd_err.pdf")
plt.clf()