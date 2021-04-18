import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# a)
data_a = pd.read_csv('build/ex2_a.csv')
ax = data_a.plot(x='x',y='rel_err_rounded',kind='scatter',marker='.',label='rel_err_rounded')
data_a.plot(x='x',y='rel_err_better',kind='scatter',marker='.',color='r',label='rel_err_better',ax=ax)
plt.savefig('build/plot_ex2_a.pdf')

# b)
plt.clf()
data_b = pd.read_csv('build/ex2_b.csv')
ax = data_b.plot(x='x',y='rel_err_rounded',kind='scatter',marker='.',label='rel_err_rounded')
data_b.plot(x='x',y='rel_err_better',kind='scatter',marker='.',color='r',label='rel_err_better',ax=ax)
plt.savefig('build/plot_ex2_b.pdf')

# b)
plt.clf()
data_c = pd.read_csv('build/ex2_c.csv')
ax = data_c.plot(x='x',y='rel_err_rounded',kind='scatter',marker='.',label='rel_err_rounded')
data_c.plot(x='x',y='rel_err_better',kind='scatter',marker='.',color='r',label='rel_err_better',ax=ax)
plt.savefig('build/plot_ex2_c.pdf')