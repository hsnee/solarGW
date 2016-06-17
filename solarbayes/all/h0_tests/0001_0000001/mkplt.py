import numpy as np
import matplotlib.pyplot as plt

p = np.loadtxt('p.txt')
p = np.exp(p-np.max(p))
a = np.loadtxt('h0.txt')

plt.figure()
plt.plot(a,p,'+')
plt.title('Probability distribution')
plt.xlabel('h0')
plt.ylabel('p')
plt.savefig('prob.png')

