import numpy as np 
import matplotlib.pyplot as plt


time = np.array(range(0,3001,100))
dt = []; dt0 = 0.1;
while dt0 < 15:
	dt.append(dt0)
	dt0 *= 1.5

l2_0 = []
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t0.100000.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t0.150000.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t0.225000.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t0.337500.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t0.506250.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t0.759375.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t1.139063.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t1.708594.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t2.562891.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t3.844336.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t5.766504.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t8.649756.txt'))
l2_0.append(np.loadtxt('../sols/stab/l2_a0.000000_t12.974634.txt'))


l2_1 = []
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t0.100000.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t0.150000.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t0.225000.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t0.337500.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t0.506250.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t0.759375.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t1.139063.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t1.708594.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t2.562891.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t3.844336.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t5.766504.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t8.649756.txt'))
l2_1.append(np.loadtxt('../sols/stab/l2_a0.500000_t12.974634.txt'))

l2_2 = []
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t0.100000.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t0.150000.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t0.225000.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t0.337500.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t0.506250.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t0.759375.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t1.139063.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t1.708594.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t2.562891.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t3.844336.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t5.766504.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t8.649756.txt'))
l2_2.append(np.loadtxt('../sols/stab/l2_a1.000000_t12.974634.txt'))

f = plt.figure(figsize=(12,12))

f.add_subplot(3,1,1)
for i in range(len(l2_0)):
	plt.semilogy(time, l2_0[i], label = 'dt=%2.2f s'%(dt[i]))
plt.xlabel('Iteration')
plt.ylabel('L2 norm')
plt.grid()
plt.legend(prop={'size': 9})
plt.ylim(0, 2)
plt.xlim(0,3500)
plt.title('L2 norm of transient solutions at various time step, when alpha = 0')


f.add_subplot(3,1,2)
for i in range(len(l2_1)):
	plt.semilogy(time, l2_1[i], label = 'dt=%2.2f s'%(dt[i]))
plt.xlabel('Iteration')
plt.ylabel('L2 norm')
plt.grid()
plt.legend(prop={'size': 9})
plt.ylim(0, 2)
plt.xlim(0,3500)
plt.title('L2 norm of transient solutions at various time step, when alpha = 0.5')

f.add_subplot(3,1,3)
for i in range(len(l2_2)):
	plt.semilogy(time, l2_2[i], label = 'dt=%2.2f s'%(dt[i]))
plt.xlabel('Iteration')
plt.ylabel('L2 norm')
plt.grid()
plt.legend(prop={'size': 9})
plt.ylim(0, 2)
plt.xlim(0,3500)
plt.title('L2 norm of transient solutions at various time step, when alpha = 1')


plt.subplots_adjust(wspace = .3, hspace = .3)
plt.show(block = False)
plt.savefig('../sols/stab/L2stab.png')
plt.close()


'''
f = plt.figure(figsize=(16,8))
plt.plot([0, 0.5, 1], [l2[0][-1], l2[1][-1], l2[2][-1]], marker='o')
plt.xlabel('alpha')
plt.ylabel('L2 norm')
plt.title('L2 norm of different alpha at t = 3000 s')
plt.grid()
plt.show(block = False)
plt.savefig('../sols/L2norm_t3000.png')
plt.close()
'''