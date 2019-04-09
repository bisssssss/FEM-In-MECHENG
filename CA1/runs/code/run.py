import os
import matplotlib.pyplot as plt
import numpy as np

def main():
	
	os.system("rm CMakeCache.txt cmake_install.cmake Makefile main1")
	os.system("rm -r CMakeFiles")
	
	os.system("cmake CMakeLists.txt")
	
	os.system("make run")
	
	os.system("mv *.h5 ../sols/h5/")
	
	# comparision
	os.system("diff ../sols/vtk/solution_p1_10_linear.vtk ../sols_sample/solution_i_linear.vtk")
	os.system("diff ../sols/vtk/solution_p1_10_quadratic.vtk ../sols_sample/solution_i_quadratic.vtk")
	os.system("diff ../sols/vtk/solution_p1_10_cubic.vtk ../sols_sample/solution_i_cubic.vtk")
	os.system("diff ../sols/vtk/solution_p2_10_linear.vtk ../sols_sample/solution_ii_linear.vtk")
	os.system("diff ../sols/vtk/solution_p2_10_quadratic.vtk ../sols_sample/solution_ii_quadratic.vtk")
	#os.system("diff ../sols/vtk/solution_p2_10_cubic.vtk solution_ii_cubic.vtk")
	
	plterr();

def plterr():
	err = np.loadtxt('../sols/plots/data/l2err.txt')
	err.reshape(4, 12)
	nmesh = np.array((3,10,100,1000)) / 0.1
	f = plt.figure(figsize = (16,12))
	for i in range(12):
		f.add_subplot(4, 3, i+1)
		yv = np.zeros(4)
		for j in range(4):
			yv[j] = err[j][i]
		slope = (np.log(yv[2]) - np.log(yv[1])) / (np.log(nmesh[2]) - np.log(nmesh[1]))
		plt.loglog(nmesh, yv, linewidth = 2, color = 'blue')
		plt.legend(['slope = %.2f'%(slope)])
		plt.title(plttitle(i) + 'L2error')
		plt.xlabel(r'1/$h_e$')
		plt.ylabel(r'error')
		plt.grid()

	plt.subplots_adjust(wspace = .3, hspace = .5)
	f.patch.set_facecolor('white')
	plt.show( block = False )
	plt.close(f)
	f.savefig('../sols/plots/L2error.png')

	err = np.loadtxt('../sols/plots/data/h1err.txt')
	err.reshape(4, 12)

	f = plt.figure(figsize = (16,12))
	for i in range(12):
		f.add_subplot(4, 3, i+1)
		yv = np.zeros(4)
		for j in range(4):
			yv[j] = err[j][i]
		slope = (np.log(yv[3]) - np.log(yv[1])) / (np.log(nmesh[3]) - np.log(nmesh[1]))
		plt.loglog(nmesh, yv, linewidth = 2, color = 'blue')
		plt.legend(['slope = %.2f'%(slope)])
		plt.title(plttitle(i) + 'H1error')
		plt.xlabel(r'1/$h_e$')
		plt.ylabel(r'error')
		plt.grid()

	plt.subplots_adjust(wspace = .3, hspace = .5)
	f.patch.set_facecolor('white')
	plt.show( block = False )
	plt.close(f)
	f.savefig('../sols/plots/H1error.png')


def plttitle(i):
	a = i // 3
	b = i % 3
	if a == 0 :
		fname = 'prob1 '
	elif a == 1:
		fname = 'prob2 '
	elif a == 2:
		fname = 'prob3 '
	elif a == 3:
		fname = 'prob4 '
	if b == 0:
		fname += 'linear '
	elif b == 1:
		fname += 'quadratic '
	elif b == 2:
		fname += 'cubic '

	return fname

if __name__ == "__main__":
	main()
