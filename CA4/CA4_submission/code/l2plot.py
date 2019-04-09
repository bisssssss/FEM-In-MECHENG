import matplotlib.pyplot as plt
import numpy as np 


l2 = []
l2.append(np.loadtxt('../sols/l2_0.txt'))
l2.append(np.loadtxt('../sols/l2_0.5.txt'))
l2.append(np.loadtxt('../sols/l2_1.txt'))

time = np.array(range(0,3001,100))

f = plt.figure(figsize=(12,12))

f.add_subplot(3,1,1)
plt.semilogy(time, l2[0], label = 'Alpha 0')
plt.xlabel('time [sec]')
plt.ylabel('L2 norm')
plt.grid()
plt.title('L2 norm of transient solutions when alpha = 0')

f.add_subplot(3,1,2)
plt.semilogy(time, l2[1], label = 'Alpha 0.5')
plt.xlabel('time [sec]')
plt.ylabel('L2 norm')
plt.grid()
plt.title('L2 norm of transient solutions when alpha = 0.5')

f.add_subplot(3,1,3)
plt.semilogy(time, l2[2], label = 'Alpha 1')
plt.xlabel('time [sec]')
plt.ylabel('L2 norm')
plt.grid()
plt.title('L2 norm of transient solutions when alpha = 1')

plt.subplots_adjust(wspace = .3, hspace = .5)
plt.show(block = False)
plt.savefig('../sols/L2plot.png')
plt.close()



f = plt.figure(figsize=(16,8))
plt.plot([0, 0.5, 1], [l2[0][-1], l2[1][-1], l2[2][-1]], marker='o')
plt.xlabel('alpha')
plt.ylabel('L2 norm')
plt.title('L2 norm of different alpha at t = 3000 s')
plt.grid()
plt.show(block = False)
plt.savefig('../sols/L2norm_t3000.png')
plt.close()