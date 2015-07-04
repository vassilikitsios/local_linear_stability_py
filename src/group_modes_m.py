# module navier_stokes_m.py

import math
import numpy as np
import scipy as sci

import io_m as io
import os

#=======================================================================
def group_modes(Nev,N_kx_real,N_kx_imag,N_kz_real,kx,kz,num_modes_to_write):
	if not os.path.exists('./results.groups'):
        	os.mkdir('./results.groups')
	if not os.path.exists('./images.groups'):
        	os.mkdir('./images.groups')

	w_all = np.array(np.zeros((Nev,N_kx_real,N_kx_imag,N_kz_real), dtype=np.complex))
	io.read_all_eigenvalues(Nev,N_kx_real,N_kx_imag,N_kz_real,w_all)

	w_grouped = np.array(np.zeros((Nev,N_kx_real,N_kx_imag,N_kz_real), dtype=np.complex))
	N_groups = group_eigenvalues(Nev,N_kx_real,N_kx_imag,N_kz_real,kx,w_all,w_grouped)

	if (N_kz_real==1):
		num_points=100
		points = np.array(np.zeros((N_kx_real*N_kx_imag,2), dtype=np.float))
		valuesI = np.imag(kx[:,:]).flatten()
		valuesR = np.real(kx[:,:]).flatten()
		for k in range(0,N_groups):
			w_imag_min = np.min(np.min(np.imag(w_grouped[k,:,:,0])))
			w_imag_max = np.max(np.max(np.imag(w_grouped[k,:,:,0])))
			if ( (w_imag_min<0.0) and (w_imag_max>0.0) ):
				w_real_min = np.min(np.min(np.real(w_grouped[k,:,:,0])))
				w_real_max = np.max(np.max(np.real(w_grouped[k,:,:,0])))
				frequency, zero_temporal_growth = np.mgrid[w_real_min:w_real_max:num_points*1j, 0:0:1j]
				points[:,0] = np.real(w_grouped[k,:,:,0]).flatten()
				points[:,1] = np.imag(w_grouped[k,:,:,0]).flatten()
				growth_rate = griddata(points, valuesI, (frequency, zero_temporal_growth), method='cubic')
				wavenumber  = griddata(points, valuesR, (frequency, zero_temporal_growth), method='cubic')
				io.write_spatial_growth_rates(frequency[1:-2],growth_rate[1:-2],wavenumber[1:-2],\
					'./results.groups/spatial_growth.mode_group'+'%04d'%(k+1)+'.dat')
			        io.plot_profile(-growth_rate,frequency/2.0/math.pi,'f/2/pi,','-imag(kx)',\
					'./images.groups/spatial_growth.mode_group'+'%04d'%(k+1)+'.kxI.png')
			        io.plot_profile(wavenumber,frequency/2.0/math.pi,'f/2/pi,','real(kx)',\
					'./images.groups/spatial_growth.mode_group'+'%04d'%(k+1)+'.kxR.png')

	io.write_grouped_modes(N_groups,N_kx_real,N_kx_imag,N_kz_real,kx,kz,w_grouped,num_modes_to_write)
	io.plot_grouped_modes(N_groups,N_kx_real,N_kx_imag,N_kz_real,kx,kz,w_grouped,num_modes_to_write)
	for k in range(0,N_groups):
                prefix = '.mode_group'+'%04d'%(k+1)
		if (N_kz_real==1):
	        	io.plot_eigenvalues(w_grouped[k,:,:,0],'./images.groups/eigenvalues'+prefix+'.png')
        		io.plot_eigenvalues(w_grouped[k,:,:,0]/np.real(kx[:,:]),'./images.groups/wave_speed'+prefix+'.png')
		elif (N_kx_imag==1):
	        	io.plot_eigenvalues(w_grouped[k,:,0,:],'./images.groups/eigenvalues'+prefix+'.png')
        		io.plot_eigenvalues(w_grouped[k,:,0,:]/np.real(kx[:,0]),'./images.groups/wave_speed'+prefix+'.png')
		elif (N_kx_real==1):
	        	io.plot_eigenvalues(w_grouped[k,0,:,:],'./images.groups/eigenvalues'+prefix+'.png')
        		io.plot_eigenvalues(w_grouped[k,0,:,:]/np.real(kx[0,:]),'./images.groups/wave_speed'+prefix+'.png')
	return

#=======================================================================
def group_eigenvalues(Nev,N_kx_real,N_kx_imag,N_kz_real,kx,w_all,w_grouped):
	print "   Grouping eigenvalues "

	BIG_NUM = 1.0e9
	SMALL_NUM = 1.0e-9

	# First set zero eigenvalues and continuous spectrum to large number
        for k in range(0,Nev):
		for i in range(0,N_kx_real):
                	for j in range(0,N_kx_imag):
                       		for m in range(0,N_kz_real):
					if ( (np.real(w_all[k,i,j,m]/kx[i,j])<=SMALL_NUM) or (np.real(w_all[k,i,j,m]/kx[i,j])>=0.9) ):
						w_all[k,i,j,m] = BIG_NUM

	# Only group discrete modes, that is those not part of the continuous spectrum
	c_diff = np.array(np.zeros((Nev), dtype=np.complex))
	counter = 0
        for k in range(0,Nev):
		if ( (np.real(w_all[k,0,0,0]/kx[0,0])<0.9) and (np.real(w_all[k,0,0,0]/kx[0,0])>SMALL_NUM) and (w_all[k,0,0,0]<BIG_NUM) ):
			w_grouped[counter,0,0,0] = w_all[k,0,0,0]
		        for i in range(0,N_kx_real):
                		for j in range(0,N_kx_imag):
                        		for m in range(0,N_kz_real):
						if ( (i>0) or (j>0) or (m>0) ):
							if (i>0):
								i0=i-1
								j0=j
								m0=m
							else:
								i0=i
								j0=max(0,j-1)
								m0=max(0,m-1)
							#c_diff[:] = w_grouped[counter,i0,j0,m0]/np.real(kx[i0,j0])*np.ones(Nev) \
							#		- w_all[:,i,j,m]/np.real(kx[i,j])
							c_diff[:] = w_grouped[counter,i0,j0,m0]/kx[i0,j0]*np.ones(Nev)-w_all[:,i,j,m]/kx[i,j]
							c_diff[:] = c_diff[:] * c_diff[:].conj()
							ev_num = np.argmin(np.real(c_diff))
							w_grouped[counter,i,j,m] = w_all[ev_num,i,j,m]
							w_all[ev_num,i,j,m] = BIG_NUM
			counter = counter + 1
	print "      number of mode groups = ", counter
	return counter

#=======================================================================
