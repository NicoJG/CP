import numpy as np
import matplotlib.pyplot as plt

rel_err_norm_a,rel_err_sym_a,rel_err_norm_b,rel_err_sym_b = np.genfromtxt("build/ex3_results.csv",delimiter=',',unpack=True)

n = np.arange(0,len(rel_err_norm_a))

plt.subplots(1,2)

plt.subplot(1,2,1)
plt.plot(n,rel_err_norm_a,'.',label='rel_err_norm_a')
plt.plot(n,rel_err_sym_a,'.',label='rel_err_sym_a')
plt.legend()

plt.subplot(1,2,2)
plt.plot(n,rel_err_norm_b,'.',label='rel_err_norm_b')
plt.plot(n,rel_err_sym_b,'.',label='rel_err_sym_b')

plt.legend()
plt.show()