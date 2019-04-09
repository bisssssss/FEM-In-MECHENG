import matplotlib.pyplot as plt
import numpy as np 


l2 = []
l2.append(np.loadtxt('../sols/l2_0.txt'))
l2.append(np.loadtxt('../sols/l2_0.5.txt'))
l2.append(np.loadtxt('../sols/l2_1.txt'))

l2lumped = []
l2lumped.append(np.loadtxt('../sols/lumped/l2_lumped_0.txt'))

l2lumped.append(np.loadtxt('../sols/lumped/l2_lumped_0.5.txt'))
l2lumped.append(np.loadtxt('../sols/lumped/l2_lumped_1.txt'))


time = np.array(range(0,3001,100))

f = plt.figure(figsize=(12,12))

f.add_subplot(3,1,1)
plt.semilogy(time, l2[0], label = 'consistent M', marker='*')
plt.semilogy(time, l2lumped[0], label = 'lumped M', marker='o')
plt.xlabel('time [sec]')
plt.ylabel('L2 norm')
plt.legend()
plt.grid()
plt.title('L2 norm of transient solutions when alpha = 0')

f.add_subplot(3,1,2)
plt.semilogy(time, l2[1], label = 'consistent M', marker='*')
plt.semilogy(time, l2lumped[1], label = 'lumped M', marker='o')
plt.xlabel('time [sec]')
plt.ylabel('L2 norm')
plt.legend()
plt.grid()
plt.title('L2 norm of transient solutions when alpha = 0.5')

f.add_subplot(3,1,3)
plt.semilogy(time, l2[2], label = 'consistent M', marker='*')
plt.semilogy(time, l2lumped[2], label = 'lumped M', marker='o')
plt.xlabel('time [sec]')
plt.ylabel('L2 norm')
plt.legend()
plt.grid()
plt.title('L2 norm of transient solutions when alpha = 1')

plt.subplots_adjust(wspace = .3, hspace = .5)
plt.show(block = False)
plt.savefig('../sols/lumped/L2plot.png')
plt.close()



f = plt.figure(figsize=(16,8))
plt.plot([0, 0.5, 1], [l2[0][-1], l2[1][-1], l2[2][-1]], label = 'consistent M', marker='*')
plt.plot([0, 0.5, 1], [l2lumped[0][-1], l2lumped[1][-1], l2lumped[2][-1]], label = 'lumped M', marker='o')
plt.xlabel('alpha')
plt.ylabel('L2 norm')
plt.title('L2 norm of different alpha at t = 3000 s')
plt.grid()
plt.legend()
plt.show(block = False)
plt.savefig('../sols/lumped/L2norm_t3000.png')
plt.close()
