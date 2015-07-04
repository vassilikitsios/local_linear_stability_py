# module io_m.py

import math
import string
import sys
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

#=======================================================================
def write_base_flow(eta,u,u_y,u_yy):
        filename = "./results/base_flow.dat"
        print '\n   Writing base flow to ', filename
        file = open(filename,'w')
        file.write('#\n')
        file.write('# y, u, u_y, u_yy\n')
        file.write('#\n')
        for i in range(0,len(eta)):
                file.write('%18.10e %18.10e %18.10e %18.10e\n' % (eta[i], u[i], u_y[i], u_yy[i]) )
        file.close()
        return

#=======================================================================
def write_eigenvalues(ev,kx,kz,N,filename):
	file = open(filename,'w')
	file.write('#\n')
	file.write('# ev num, ev_r, ev_i, c_r, c_i, kx_r, kx_i, kz_r, kz_i\n')
	file.write('#\n')
	for j in range(0,N):
        	file.write('%4.1d %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e\n' % \
			(j+1, ev[j].real, ev[j].imag, np.real(ev[j]/kx), np.imag(ev[j]/kx), kx.real, kx.imag, kz.real, kz.imag) )
	file.close()
	return

#=======================================================================
def write_eigenvectors(eta,ev,vr,vl,kx,kz,N_modes,prefix):
	N = len(eta)-1
	if (len(vr)>2*(N+1)):
		write_eigenvector(eta,ev,vr[0      :   N+1 ],vl[0      :   N+1 ],kx,kz,N_modes,prefix+'.u.dat')
		write_eigenvector(eta,ev,vr[N+1    :2*(N+1)],vl[N+1    :2*(N+1)],kx,kz,N_modes,prefix+'.v.dat')
		write_eigenvector(eta,ev,vr[2*(N+1):3*(N+1)],vl[2*(N+1):3*(N+1)],kx,kz,N_modes,prefix+'.w.dat')
		write_eigenvector(eta,ev,vr[3*(N+1):4*(N+1)],vl[3*(N+1):4*(N+1)],kx,kz,N_modes,prefix+'.p.dat')
	else:
		write_eigenvector(eta,ev,vr[0      :   N+1 ],vl[0      :   N+1 ],kx,kz,N_modes,prefix+'.v.dat')
		write_eigenvector(eta,ev,vr[N+1    :2*(N+1)],vl[N+1    :2*(N+1)],kx,kz,N_modes,prefix+'.vor.dat')
	return

#=======================================================================
def write_eigenvector(eta,ev,vr,vl,kx,kz,N_modes,prefix):
	N = len(eta)-1
	for j in range(0,N_modes):
		filename = prefix + '.mode' + '%04d'%(j+1) + '.dat'
		file = open(filename,'w')
		file.write('#\n')
		file.write('# kx = %18.10e %18.10e\n' % (kx.real, kx.imag) )
		file.write('# kz = %18.10e %18.10e\n' % (kz.real, kz.imag) )
		file.write('# evalue = %18.10e %18.10e\n' % (ev[j].real, ev[j].imag) )
		file.write('# wave speed = %18.10e %18.10e\n' % (np.real(ev[j]/kx), np.imag(ev[j]/kx)) )
		file.write('#\n')
		file.write('# eta, ev_r, ev_i, evA_r, evA_i\n')
		file.write('#\n')
		for i in range(0,N+1):
        		file.write('%18.10e %18.10e %18.10e %18.10e %18.10e\n' % (eta[i],vr[i,j].real,vr[i,j].imag,vl[i,j].real,vl[i,j].imag))
        		#file.write('%18.10e %18.10e %18.10e %18.10e %18.10e\n' % (eta[i],vr[j,i].real,vr[j,i].imag,vl[j,i].real,vl[j,i].imag))
		file.close()
	return

#=======================================================================
def read_all_eigenvalues(Nev,N_kx_real,N_kx_imag,N_kz_real,w_all):
        w = np.array(np.zeros((Nev), dtype=np.complex))
        print "   Reading eigenvalues"
        for i in range(0,N_kx_real):
                for j in range(0,N_kx_imag):
                        for m in range(0,N_kz_real):
                                filename = './results/eigenvalues.kxR'+'%04d'%(i+1)+'.kxI'+'%04d'%(j+1)+'.kzR'+'%04d'%(m+1)+'.dat'
                                read_eigenvalues(Nev,filename,w)
                                w_all[:,i,j,m] = w[:]
        return

#=======================================================================
def read_eigenvalues(Nev,filename,w):
        file = open(filename,'r')
        for k in range(0,3):
                linestring = file.readline()
        for k in range(0,Nev):
                linestring = file.readline()
                linelist = string.split(linestring)
                w_real = float(linelist[1])
                w_imag = float(linelist[2])
                w[k] = w_real + 1j*w_imag
        file.close()
        return

#=======================================================================
def write_grouped_modes(N_groups,N_kx_real,N_kx_imag,N_kz_real,kx,kz,w_grouped,num_modes_to_write):
        print "   Writing grouped eigenvalues"
        for k in range(0,min(num_modes_to_write,N_groups)):
                filename = './results.groups/eigenvalues.mode' + '%04d'%(k+1) + '.dat'
                file = open(filename,'w')
                file.write('#\n')
                file.write('# kx_r-index, kx_i-index, kz_r-index, ev_r, ev_i, c_r, c_i, kx_r, kx_i, kz_r, kz_i\n')
                file.write('#\n')
                for i in range(0,N_kx_real):
                        for j in range(0,N_kx_imag):
                                for m in range(0,N_kz_real):
                                        file.write('%4.1d %4.1d %4.1d %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e %18.10e\n' % \
                                                (i,j,m,w_grouped[k,i,j,m].real, w_grouped[k,i,j,m].imag,
						np.real(w_grouped[k,i,j,m]/kx[i,j]), np.imag(w_grouped[k,i,j,m]/kx[i,j]), \
                                                kx[i,j].real, kx[i,j].imag, kz[m].real, kz[m].imag) )
        	file.close()
        return

#=======================================================================
def write_spatial_growth_rates(frequency,growth_rate,wavenumber,filename):
	file = open(filename,'w')
	file.write('#\n')
	file.write('# ev_r, kx_i, kx_r\n')
	file.write('#\n')
	for k in range(0,len(frequency)):
		file.write('%18.10e %18.10e %18.10e\n' % (frequency[k], growth_rate[k], wavenumber[k]) )
       	file.close()

#=======================================================================
def plot_eigenvectors(N_modes,eta,vr,vl,prefix):
	for j in range(0,N_modes):
		plot_eigenvector(eta,vr[:,j],"eigenvector","eta",prefix + '.mode' + '%04d'%(j+1) + '.direct')
		plot_eigenvector(eta,vl[:,j],"eigenvector","eta",prefix + '.mode' + '%04d'%(j+1) + '.adjoint')
	return

#=======================================================================
def plot_eigenvector(y,u,xlabel,ylabel,prefix):
	N = len(y)-1
	if (len(u)>2*(N+1)):
		plot_complex_profile(y,u[0      :   N+1 ],"eigenvector","y",prefix+'.u.png')
		plot_complex_profile(y,u[N+1    :2*(N+1)],"eigenvector","y",prefix+'.v.png')
		plot_complex_profile(y,u[2*(N+1):3*(N+1)],"eigenvector","y",prefix+'.w.png')
		plot_complex_profile(y,u[3*(N+1):4*(N+1)],"eigenvector","y",prefix+'.p.png')
	else:
		plot_complex_profile(y,u[0      :   N+1 ],"eigenvector","y",prefix+'.v.png')
		plot_complex_profile(y,u[N+1    :2*(N+1)],"eigenvector","y",prefix+'.vort.png')
	return

#=======================================================================
def plot_eigenvalues(ev,filename):
	print "      ", filename
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        plt.plot(ev[:].real,ev[:].imag,'or',[-0.1,1.1],[0.0,0.0],'-k')
	plt.xlim(-0.1,1.1)
	plt.ylim(-1.1,1.1)
        #plt.legend()
        ax1.set_xlabel('real(ev)')
        ax1.set_ylabel('imag(ev)')
        fig.savefig(filename)
	return

#=======================================================================
def plot_base_flow(y,u,u_y,u_yy):
	plot_profile(y,u,'y,','u','./images/u.png')
	plot_profile(y,u_y,'y,','u_y','./images/u_y.png')
	plot_profile(y,u_yy,'y,','u_yy','./images/u_yy.png')
	return

#=======================================================================
def plot_profile(y,u,xlabel,ylabel,filename):
	print "      ", filename
	plt.close('all')
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        plt.plot(u,y,'-xr')
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)
        fig.savefig(filename)
	return

#=======================================================================
def plot_complex_profile(y,u,xlabel,ylabel,filename):
	print "      ", filename
	plt.close('all')
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        plt.plot(u.real,y,'-r',u.imag,y,'-b')
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)
        fig.savefig(filename)
	return

#=======================================================================
def plot_grouped_modes(N_groups,N_kx_real,N_kx_imag,N_kz_real,kx,kz,w_grouped,num_modes_to_write):
        print "   Plotting grouped eigenvalues"

	Nx = N_kx_real
	x_label = "$k_{x,r}$"
	if (N_kx_imag>1):
		Ny = N_kx_imag
		y_label = "$k_{x,i}$"
	elif (N_kz_real>1):
		Ny = N_kz_real
		y_label = "$k_z$"

	w_plot = np.array(np.zeros((Nx,Ny), dtype=np.complex))
        #for k in range(0,min(N_groups,num_modes_to_write)):
        for k in range(0,N_groups):
                filename_real = './images.groups/eigenvalues_real.mode' + '%04d'%(k+1) + '.ps'
                filename_imag = './images.groups/eigenvalues_imag.mode' + '%04d'%(k+1) + '.ps'
		if (N_kx_imag>1):
			print "      ", filename_real
        	        w_plot[:,:] = np.real(w_grouped[k,:,:,0]) / np.real(kx[:,:])
			plot_surface(kx[:,0].real,kx[0,:].imag,w_plot.T,x_label,y_label,filename_real)
			print "      ", filename_imag
        	        w_plot[:,:] = np.imag(w_grouped[k,:,:,0]) / np.real(kx[:,:])
			plot_surface(kx[:,0].real,kx[0,:].imag,w_plot.T,x_label,y_label,filename_imag)
		elif (N_kz_real>1):
			print "      ", filename_real
        	        w_plot[:,:] = w_grouped[k,:,0,:].real
			plot_surface(kx[:,0].real,kz[0,:].real,w_plot.T,x_label,y_label,filename_real)
			print "      ", filename_imag
        	        w_plot[:,:] = w_grouped[k,:,0,:].imag
			plot_surface(kx[:,0].real,kz[0,:].real,w_plot.T,x_label,y_label,filename_imag)
	return

#=======================================================================
def plot_surface(x,y,z,x_label,y_label,filename):
        #matplotlib.rc('font', size=32, family='sans-serif')
       	matplotlib.rc('font', size=24, family='latex')
        matplotlib.rc('text', usetex=True)
       	matplotlib.rcParams['xtick.major.pad']=12
        matplotlib.rcParams['ytick.major.pad']=12

       	plt.close('all')
        fig = plt.figure()
       	ax = fig.add_subplot(111)
        #plt.xlim(0,Nx)
        #plt.ylim(0,Ny)
       	ax.set_xlabel(x_label)
       	ax.set_ylabel(y_label)

        #z = z.clip(min_contour_val,max_contour_val)
        #lev = np.linspace(min_contour_val,max_contour_val,num_contours)
        CS = plt.contour(x,y,z,linewidths=0.5,colors='k')
        CS = plt.contour(x,y,z,[0],linewidths=4.0,colors='k')
        #CS = plt.contourf(x,y,z,cmap=plt.cm.gray_r)
        CS = plt.contourf(x,y,z)

        plt.subplots_adjust(left=0.2, bottom=0.2, right=0.85, top=0.9)
        plt.colorbar(format='%2.2g')
        #plt.colorbar(format='%2.1f') 

        #plt.show()
        fig.savefig(filename)
	return

#=======================================================================
