import numpy as np
import matplotlib.pyplot as plt

n,rel_err_norm_a,rel_err_sym_a,rel_err_norm_b,rel_err_sym_b = np.genfromtxt("build/ex3.csv",delimiter=',',unpack=True)

plt.subplots(1,2)

plt.subplot(1,2,1)
plt.xlabel('n')
plt.ylabel('relative error')
plt.plot(n[500:],rel_err_norm_a[500:],'.',label='rel_err_norm_a')
plt.plot(n[500:],rel_err_sym_a[500:],'.',label='rel_err_sym_a')
plt.legend()

plt.subplot(1,2,2)
plt.xlabel('n')
plt.ylabel('relative error')
plt.plot(n[500:],rel_err_norm_b[500:],'.',label='rel_err_norm_b')
plt.plot(n[500:],rel_err_sym_b[500:],'.',label='rel_err_sym_b')

plt.tight_layout()
plt.legend()
plt.savefig('build/plot_ex3.pdf')