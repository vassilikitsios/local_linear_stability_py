# module navier_stokes_m.py

import math
import numpy as np
import scipy as sci
import scipy.linalg as linalg

import io_m as io
import os

#=======================================================================
def generate_wavenumber_sets(kx_real_min, kx_real_max, N_kx_real, kx_imag_min, kx_imag_max, N_kx_imag,\
			        N_kz_real, kz_real_max, kz_real_min, kz_imag, kx, kz):
	d_kx_real = (kx_real_max-kx_real_min)/max(N_kx_real-1,1)
	d_kx_imag = (kx_imag_max-kx_imag_min)/max(N_kx_imag-1,1)
	d_kz_real = (kz_real_max-kz_real_min)/max(N_kz_real-1,1)
	for i in range(0,N_kx_real):
        	for j in range(0,N_kx_imag):
			kx[i,j] = (kx_real_min + d_kx_real*i) + (kx_imag_min + d_kx_imag*j)*1j
	for m in range(0,N_kz_real):
		kz[m] = (kz_real_min + d_kz_real*m) + kz_imag*1j
	return

#=======================================================================
def run_stability_analysis(Nev,N_kx_real,N_kx_imag,N_kz_real,kx,kz,equations,eta,u,u_y,u_yy,Re,D1,D2,D4,top_BC,bottom_BC,num_modes_to_write):
        A = np.array(np.zeros((Nev,Nev), dtype=np.complex))                     # linear operators
        B = np.array(np.zeros((Nev,Nev), dtype=np.complex))
        w = np.array(np.zeros((Nev), dtype=np.complex))                         # eigenvalues
        vl = np.array(np.zeros((Nev,Nev), dtype=np.complex))                    # left eigenvalues
        vr = np.array(np.zeros((Nev,Nev), dtype=np.complex))                    # right eigenvalues

        counter=0
        for i in range(0,N_kx_real):
                for j in range(0,N_kx_imag):
                        for m in range(0,N_kz_real):
                                counter = counter + 1
                                print '   case ', counter, ' of ', N_kx_real*N_kx_imag*N_kz_real, ' kx=', kx[i,j], ' kz=', kz[m]

                                build_perturbation_linear_operators(equations,u,u_y,u_yy,Re,kx[i,j],kz[m],D1,D2,D4,\
                                                                        top_BC,bottom_BC,A,B)
                                N_modes = calculate_and_sort_eigenvalues(Nev,A,B,w,vl,vr)
                                print "number of valid modes = ", N_modes, " of ", Nev

                                prefix = '.kxR'+'%04d'%(i+1)+'.kxI'+'%04d'%(j+1)+'.kzR'+'%04d'%(m+1)

                                io.write_eigenvalues(w,kx[i,j],kz[m],Nev,'./results/eigenvalues'+prefix+'.dat')
                                io.plot_eigenvalues(w,'./images/eigenvalues'+prefix+'.png')
                                io.plot_eigenvalues(w/kx[i,j],'./images/wave_speed'+prefix+'.png')

				io.write_eigenvectors(eta,w,vr,vl,kx[i,j],kz[m],min(N_modes,num_modes_to_write),\
                                                        './results/eigenvectors'+prefix)
				io.plot_eigenvectors(min(N_modes,num_modes_to_write),eta,vr,vl,'./images/eigenvectors'+prefix)
	return

#=======================================================================
def run_transient_growth_analysis(Nev,N_kx_real,N_kx_imag,N_kz_real,kx,kz,equations,eta,u,u_y,u_yy,Re,\
					D1,D2,D4,IW,top_BC,bottom_BC,num_modes_to_write,dt,num_time_steps):
        A  = np.array(np.zeros((Nev,Nev), dtype=np.complex))                     # linear operators
        B  = np.array(np.zeros((Nev,Nev), dtype=np.complex))
        w  = np.array(np.zeros((Nev), dtype=np.complex))                         # eigenvalues
        vl = np.array(np.zeros((Nev,Nev), dtype=np.complex))                     # left eigenvalues
        vr = np.array(np.zeros((Nev,Nev), dtype=np.complex))                     # right eigenvalues

	N = len(IW)-1
        IW_all = np.array(np.zeros((N+1)*4, dtype=np.float))
	IW_all[0      :N+1    ] = IW
	IW_all[N+1    :2*(N+1)] = IW
	IW_all[2*(N+1):3*(N+1)] = IW
	#IW_all[0      :N+1    ] = 1.0
	#IW_all[N+1    :2*(N+1)] = 1.0
	#IW_all[2*(N+1):3*(N+1)] = 1.0
	IW_all[3*(N+1):4*(N+1)] = 0.0

	time = np.array(np.zeros((num_time_steps), dtype=np.float))		# time series
        for t in range(0,num_time_steps):
		time[t] = dt*t

        counter=0
        for i in range(0,N_kx_real):
                for j in range(0,N_kx_imag):
                        for m in range(0,N_kz_real):
                                counter = counter + 1
                                print '   case ', counter, ' of ', N_kx_real*N_kx_imag*N_kz_real, ' kx=', kx[i,j], ' kz=', kz[m]

                                build_perturbation_linear_operators(equations,u,u_y,u_yy,Re,kx[i,j],kz[m],D1,D2,D4,\
                                                                        top_BC,bottom_BC,A,B)
                                N_modes = calculate_and_sort_eigenvalues(Nev,A,B,w,vl,vr)
                                print "number of valid modes = ", N_modes, " of ", Nev

                                prefix = '.kxR'+'%04d'%(i+1)+'.kxI'+'%04d'%(j+1)+'.kzR'+'%04d'%(m+1)
                                io.write_eigenvalues(w,kx[i,j],kz[m],Nev,'./results/eigenvalues'+prefix+'.dat')
                                io.write_eigenvectors(eta,w,vr,vl,kx[i,j],kz[m],min(N_modes,num_modes_to_write),\
                                                        './results/eigenvectors'+prefix)
                                io.plot_eigenvalues(w,'./images/eigenvalues'+prefix+'.png')
                                io.plot_eigenvalues(w/kx[i,j],'./images/wave_speed'+prefix+'.png')
				io.plot_eigenvectors(min(N_modes,num_modes_to_write),eta,vr,vl,'./images/eigenvectors'+prefix)

				#N_modes_check = N_modes
				N_modes_check = 30
				calculate_transient_growth(w[0:N_modes_check],vr[:,0:N_modes_check],IW_all,N,N_modes_check,time)
	return

#=======================================================================
def calculate_transient_growth(w,vr,IW_all,N,N_modes,time):
	
# first normalise the eigenvectors with respect to energy
# then compute the inner product of the eigenvectors in energy norm
        temp = np.array(np.zeros((4*(N+1),N_modes), dtype=np.complex))
	for i in range(0,N_modes):
		temp[:,i] = vr[:,i] * IW_all[:]
		mag = np.dot(temp[:,i].T.conj(),temp[:,i])
		temp[:,i] = temp[:,i] / mag
        inner_product = np.array(np.zeros((N_modes,N_modes), dtype=np.complex))
	inner_product = np.dot( temp.T.conj(), temp )

# compute decomposition inner_product=F^*F
        U = np.array(np.zeros((N_modes,N_modes), dtype=np.complex))
        S = np.array(np.zeros((N_modes), dtype=np.complex))
        V = np.array(np.zeros((N_modes,N_modes), dtype=np.complex))
	U, S, V = linalg.svd(inner_product, full_matrices=True)

        S_sqrt = np.array(np.zeros((N_modes), dtype=np.float))
	for i in range(0,N_modes):
		S_sqrt[i] = np.sqrt(S[i].real)

        F = np.array(np.zeros((N_modes,N_modes), dtype=np.complex))
	F = np.dot(np.diag(S_sqrt,0),U.T.conj())

        invF = np.array(np.zeros((N_modes,N_modes), dtype=np.complex))
	invF = np.dot(U,np.diag(np.ones(len(S_sqrt))/S_sqrt,0))

# compute transient growth
        E = np.array(np.zeros((N_modes,N_modes), dtype=np.complex))
        Q = np.array(np.zeros((N_modes,N_modes), dtype=np.complex))
	max_growth = 0
	max_growth_time = 0
	for t in range(0,len(time)):
		for i in range(0,N_modes):
			Q[i,i] = math.exp(-1j*w[i]*time[t])
		E = np.dot(np.dot(F,Q),invF)
		Ue, Se, Ve = np.linalg.svd(E, full_matrices=True)
		growth = math.pow(Se[0],2.0)
		print 'time, energy growth: ', time[t], growth
		if (growth>max_growth):
			max_growth = growth
			max_growth_time = time[t] 
	print 'max_growth, max_growth_time: ', max_growth, max_growth_time

	return

#=======================================================================
def calculate_and_sort_eigenvalues(Nev,A,B,w_sorted,vl_sorted,vr_sorted):
        w = np.array(np.zeros((Nev), dtype=np.complex))                         # eigenvalues
        vl = np.array(np.zeros((Nev,Nev), dtype=np.complex))                    # left eigenvalues
        vr = np.array(np.zeros((Nev,Nev), dtype=np.complex))                    # right eigenvalues
        w_imag_sorted = np.array(np.zeros((Nev), dtype=np.float))		# sort on basis of imaginary component of the eigenvalues 
	mode_index = np.zeros((Nev), dtype=np.integer)                          # sorted mode index
        
	# Calculate eigenvalues
        [w,vl,vr] = linalg.eig(A, B, left=True, right=True, overwrite_a=False, overwrite_b=False)

	# Find out the number of valid modes 
	BIG_NUMBER = 1.0e10
        N_modes = 0
        for k in range(0,Nev):
                mode_index[k] = k
                if (math.isnan(w[k].real) or math.isnan(w[k].imag) or (w[k].real<=0.0) ):
                        w_imag_sorted[k] = -BIG_NUMBER
                        w[k] = -BIG_NUMBER -BIG_NUMBER*1j
                        vl[:,k] = 0.0
                        vr[:,k] = 0.0
                else:
                        w_imag_sorted[k] = w[k].imag
                        N_modes += 1

	# Sort on basis of imaginary component of the eigenvalues
        w_imag_sorted, mode_index = zip(*sorted(zip(w_imag_sorted, mode_index),reverse=True))

	# Rearrange eigenvalues and eigenvectors
        w_sorted[:] = -BIG_NUMBER
        for k in range(0,N_modes):
                w_sorted[k] = w[mode_index[k]]
                vl_sorted[:,k] = vl[:,mode_index[k]]
                vr_sorted[:,k] = vr[:,mode_index[k]]

	return N_modes

#=======================================================================
def build_perturbation_linear_operators(equations,u,u_y,u_yy,Re,kx,kz,D1,D2,D4,top_BC,bottom_BC,A,B):
        if (equations=="velocity_formulation"):
                build_velocity_perturbation_linear_operators(u,u_y,Re,kx,kz,D1,D2,top_BC,bottom_BC,A,B)
        elif (equations=="vorticity_formulation"):
                build_vorticity_perturbation_linear_operators(u,u_y,u_yy,Re,kx,kz,D1,D2,D4,top_BC,bottom_BC,A,B)
        return

#=======================================================================
def build_vorticity_perturbation_linear_operators(u,u_y,u_yy,Re,kx,kz,D1,D2,D4,top_BC,bottom_BC,A,B):
	N = len(u)-1
	I = np.array(np.eye((N+1), dtype=np.float))
	k2 = np.real( kx*np.conj(kx) + kz*np.conj(kz) )
	build_vorticity_perturbation_linear_operator_A(I,u,u_y,u_yy,Re,kx,kz,k2,D2,D4,A)
	build_vorticity_perturbation_linear_operator_B(D2,k2,I,B)

	ev_shift = -400.0*1j
	if (bottom_BC=="wall"):	
		apply_vorticity_boundary_conditions_to_bottom_wall(D1,ev_shift,kx,kz,A,B)
	elif (bottom_BC=="freestream"):	
		apply_vorticity_boundary_conditions_to_bottom_freestream(kx,kz,D1,D2,A,B)
	if (top_BC=="wall"):	
		apply_vorticity_boundary_conditions_to_top_wall(D1,ev_shift,kx,kz,A,B)
	elif (top_BC=="freestream"):	
		apply_vorticity_boundary_conditions_to_top_freestream(kx,kz,D1,D2,A,B)
	return 

#=======================================================================
def build_vorticity_perturbation_linear_operator_B(D2,k2,I,B):
	N = len(I)-1
	B[:,:] = 0.0
	B[ 0:N+1, 0:N+1 ] = -1j*(k2*I-D2)
	B[ N+1:2*(N+1), N+1:2*(N+1) ] = -1j*I
	return

#=======================================================================
def build_vorticity_perturbation_linear_operator_A(I,u,u_y,u_yy,Re,kx,kz,k2,D2,D4,A):
	N = len(u)-1
	diag_u = np.array(np.zeros((N+1,N+1), dtype=np.float))
	diag_u_y = np.array(np.zeros((N+1,N+1), dtype=np.float))
	diag_u_yy = np.array(np.zeros((N+1,N+1), dtype=np.float))
	for i in range(0,N+1):
		diag_u[i,i] = u[i]
		diag_u_y[i,i] = u_y[i]
		diag_u_yy[i,i] = u_yy[i]
	
	A[:,:] = 0.0
	A[ 0:N+1, 0:N+1 ] = -1j*kx*np.dot(diag_u,k2*I-D2) - 1j*kx*diag_u_yy - 1.0/Re*(D4 - 2.0*D2*k2 + k2*k2*I)
	A[ 0:N+1, N+1:2*(N+1) ] = 0.0
	
	A[ N+1:2*(N+1), 0:N+1 ] = -1j*kz*diag_u_y
	A[ N+1:2*(N+1), N+1:2*(N+1) ] = -1j*kx*diag_u - 1.0/Re*(k2*I-D2)
	return	

#=======================================================================
def apply_vorticity_boundary_conditions_to_bottom_wall(D1,ev_shift,kx,kz,A,B):
	N = len(D1)-1

	# BC: real[v(xi=1)=0] + imag[v(xi=1)=0]
	A[0, :] = 0.0
	A[0, 0] = 1.0

	# BC: d/dn real[v(xi=1)] + imag[v(xi=1)] = 0
	A[1, :] = 0.0
	A[1, 0:N+1] = D1[0,:]
	
	# BC: real[zeta(xi=1)] + imag[zeta(xi=1)] = 0
	A[N+1, :] = 0.0
	A[N+1, N+1] = 1.0
	
	# Shift eigenvalues
	B[0, :] = A[0, :] / ev_shift
	B[1, :] = A[1, :] / ev_shift
	B[N+1, :] = A[N+1, :] / ev_shift

	return

#============================================================================
def apply_vorticity_boundary_conditions_to_top_wall(D1,ev_shift,kx,kz,A,B):
	N = len(D1)-1

	# BC: real[v(xi=-1)] + imag[v(xi=-1)] = 0
	A[N, :] = 0.0
	A[N, N] = 1.0
	
	# BC: d/dn real[v(xi=-1)] + imag[v(xi=-1)] = 0
	A[N-1, :] = 0.0
	A[N-1, 0:N+1] = D1[N,:]
	
	# BC: real[zeta(xi=-1)] + imag[zeta(xi=-1)] = 0
	A[2*(N+1)-1, :] = 0.0
	A[2*(N+1)-1, 2*(N+1)-1] = 1.0
	
	# Shift eigenvalues
	B[N, :] = A[N, :] / ev_shift
	B[N-1, :] = A[N-1, :] / ev_shift
	B[2*(N+1)-1, :] = A[2*(N+1)-1, :] / ev_shift

	return

#============================================================================
def build_velocity_perturbation_linear_operators(u,u_y,Re,kx,kz,D1,D2,top_BC,bottom_BC,A,B):
	N = len(u)-1
	I = np.array(np.eye((N+1), dtype=np.float))
	build_velocity_perturbation_linear_operator_A(I,u,u_y,Re,kx,kz,D1,D2,A)
	build_velocity_perturbation_linear_operator_B(I,B)

	ev_shift = -400.0*1j	
	if (bottom_BC=="wall"):	
		apply_velocity_boundary_conditions_to_bottom_wall(D1,ev_shift,kx,kz,A,B)
	elif (bottom_BC=="freestream"):	
		apply_velocity_boundary_conditions_to_bottom_freestream(kx,kz,D1,D2,A,B)

	if (top_BC=="wall"):	
		apply_velocity_boundary_conditions_to_top_wall(D1,ev_shift,kx,kz,A,B)
	elif (top_BC=="freestream"):	
		apply_velocity_boundary_conditions_to_top_freestream(kx,kz,D1,D2,A,B)

	return 

#=======================================================================
def build_velocity_perturbation_linear_operator_B(I,B):
	N = len(I)-1
	B[:,:] = 0.0
	B[ 0:N+1, 0:N+1 ] = -1j * I
	B[ N+1:2*(N+1), N+1:2*(N+1) ] = -1j * I
	B[ 2*(N+1):3*(N+1), 2*(N+1):3*(N+1) ] = -1j * I
	B[ 3*(N+1):4*(N+1), 3*(N+1):4*(N+1) ] = 0.0
	return

#=======================================================================
def build_velocity_perturbation_linear_operator_A(I,u,u_y,Re,kx,kz,D1,D2,A):
	N = len(D1)-1
	diag_u = np.array(np.zeros((N+1,N+1), dtype=np.float))
	diag_u_y = np.array(np.zeros((N+1,N+1), dtype=np.float))
	for i in range(0,N+1):
		diag_u[i,i] = u[i]
		diag_u_y[i,i] = u_y[i]
	k2 = kx*np.conj(kx) + kz*np.conj(kz)
	Z = np.array(np.zeros((N+1,N+1), dtype=np.complex))
	Z = -1j*kx*diag_u + 1.0/Re*(D2 - k2*I)
	
	A[:,:] = 0.0
	A[ 0:N+1, 0:N+1 ] = Z
	A[ 0:N+1, N+1:2*(N+1) ] = -diag_u_y
	A[ 0:N+1, 2*(N+1):3*(N+1) ] = 0.0
	A[ 0:N+1, 3*(N+1):4*(N+1) ] = -1j*kx * I
	
	A[ N+1:2*(N+1), 0:N+1 ] = 0.0
	A[ N+1:2*(N+1), N+1:2*(N+1) ] = Z
	A[ N+1:2*(N+1), 2*(N+1):3*(N+1) ] = 0.0
	A[ N+1:2*(N+1), 3*(N+1):4*(N+1) ] = -D1
	
	A[ 2*(N+1):3*(N+1), 0:N+1 ] = 0.0
	A[ 2*(N+1):3*(N+1), N+1:2*(N+1) ] = 0.0
	A[ 2*(N+1):3*(N+1), 2*(N+1):3*(N+1) ] = Z
	A[ 2*(N+1):3*(N+1), 3*(N+1):4*(N+1) ] = -1j*kz * I
	
	A[ 3*(N+1):4*(N+1), 0:N+1 ] = 1j*kx * I
	A[ 3*(N+1):4*(N+1), N+1:2*(N+1) ] = D1
	A[ 3*(N+1):4*(N+1), 2*(N+1):3*(N+1) ] = 1j*kz * I
	A[ 3*(N+1):4*(N+1), 3*(N+1):4*(N+1) ] = 0.0
	return

#============================================================================
def apply_velocity_boundary_conditions_to_bottom_wall(D1,ev_shift,kx,kz,A,B):
	N = len(D1)-1

	# BC: real[u(xi=1)=0] + imag[u(xi=1)=0]
	A[0, : ] = 0.0
	A[0, 0] = 1.0
	B[0, : ] = 0.0		
	
	# BC: real[v(xi=1)=0] + imag[v(xi=1)=0]
	A[N+1, : ] = 0.0
	A[N+1, N+1] = 1.0
	B[N+1, : ] = 0.0
	
	# BC: real[w(xi=1)=0] + imag[w(xi=1)=0]
	A[2*(N+1), : ] = 0.0
	A[2*(N+1), 2*(N+1)] = 1.0
	B[2*(N+1), : ] = 0.0
	
	# BC: real[dpi_dxj(xi=1)] + imag [dpi_dxj(xi=1)] = 0
	A[3*(N+1), : ] = 0.0
	A[3*(N+1), 3*(N+1):4*(N+1) ] = D1[0,:]
	#A[3*(N+1), 4*(N+1)-1 ] += 1j*kx + 1j*kz
	B[3*(N+1), : ] = 0.0
	
	# Shift eigenvalues
	B[0, : ] = A[0, : ] / ev_shift
	B[N+1, : ] = A[N+1, : ] / ev_shift
	B[2*(N+1), : ] = A[2*(N+1), : ] / ev_shift
	B[3*(N+1), : ] = A[3*(N+1), : ] / ev_shift

	return

#============================================================================
def apply_velocity_boundary_conditions_to_top_wall(D1,ev_shift,kx,kz,A,B):
	N = len(D1)-1

	# BC: real[u(xi=-1)=0] + imag[u(xi=-1)=0]
	A[N, : ] = 0.0
	A[N, N] = 1.0
	B[N, : ] = 0.0
	
	# BC: real[v(xi=-1)=0] + imag[v(xi=-1)=0]
	A[2*(N+1)-1, : ] = 0.0
	A[2*(N+1)-1, 2*(N+1)-1] = 1.0
	B[2*(N+1)-1, : ] = 0.0
	
	# BC: real[w(xi=-1)=0] + imag[w(xi=-1)=0]
	A[3*(N+1)-1, : ] = 0.0
	A[3*(N+1)-1, 3*(N+1)-1] = 1.0
	B[3*(N+1)-1, : ] = 0.0
	
	# BC: real[dpi_dxj(xi=-1)] + imag [dpi_dxj(xi=-1)] = 0
	A[4*(N+1)-1, : ] = 0.0
	A[4*(N+1)-1, 3*(N+1):4*(N+1) ] = D1[N,:]
	#A[4*(N+1)-1, 4*(N+1)-1 ] += 1j*kx + 1j*kz
	B[4*(N+1)-1, : ] = 0.0
	
	# Shift eigenvalues
	B[N, : ] = A[N, : ] / ev_shift
	B[2*(N+1)-1, : ] = A[2*(N+1)-1, : ] / ev_shift
	B[3*(N+1)-1, : ] = A[3*(N+1)-1, : ] / ev_shift
	B[4*(N+1)-1, : ] = A[4*(N+1)-1, : ] / ev_shift

	return

#============================================================================
def apply_velocity_boundary_conditions_to_bottom_freestream(kx,kz,D1,D2,A,B):
        k2 = np.real(kx*kx.conj() + kz*kz.conj())
        k = math.sqrt(k2) 
	N = len(D1)-1
        
        # BC: real[du_dy(xi=1) = -k u(xi=1)]                                            free stream
        A[0, 0:4*(N+1)-1] = 0.0
        A[0, 0:N+1] = D1[0,0:N+1]
        A[0, 0] = A[0, 0] + k
        B[0, 0:4*(N+1)-1] = 0.0

        # BC: real[dv_dy(xi=1) = -k u(xi=1)]                                            free stream
        A[N+1, 0:4*(N+1)-1] = 0.0
        A[N+1, N+1:2*(N+1)] = D1[0,0:N+1]
        A[N+1, N+1] = A[N+1, N+1] + k
        B[N+1, 0:4*(N+1)-1] = 0.0

        # BC: real[dw_dy(xi=1) = -k u(xi=1)]                                            free stream
        A[2*(N+1), 0:4*(N+1)-1] = 0.0
        A[2*(N+1), 2*(N+1):3*(N+1) ] = D1[0,0:N+1]
        A[2*(N+1), 2*(N+1)] = A[2*(N+1), 2*(N+1)] + k
        B[2*(N+1), 0:4*(N+1)-1] = 0.0

        # BC: real[d2pi_dxj2(xi=1)] = 0                                                 freestream
        A[3*(N+1), 0:4*(N+1)-1] = 0.0
        A[3*(N+1), 3*(N+1):4*(N+1) ] = D2[0,0:N+1]
        A[3*(N+1), 3*(N+1)] = A[3*(N+1), 3*(N+1)] - k2
        B[3*(N+1), 0:4*(N+1)-1] = 0.0

        return

#============================================================================
def apply_velocity_boundary_conditions_to_top_freestream(kx,kz,D1,D2,A,B):
        # Note sign of 'k' is changed to ensure the BC decay's to zero
        k2 = np.real(kx*kx.conj() + kz*kz.conj())
        k = math.sqrt(k2)
	N = len(D1)-1

        # BC: real[du_dy(xi=-1) = -k u(xi=-1)] + imag[du_dy(xi=-1) = -k u(xi=-1)]       free stream
        #A[N, 0:4*(N+1)-1] = 0.0
        A[N, :] = 0.0
        A[N, 0:N+1] = D1[N, 0:N+1]
        A[N, N] = A[N, N] + k
        B[N, :] = 0.0

        # BC: real[dv_dy(xi=-1) = -k u(xi=-1)] + imag[dv_dy(xi=-1) = -k u(xi=-1)]       free stream
        A[2*(N+1)-1, :] = 0.0
        A[2*(N+1)-1, N+1:2*(N+1) ] = D1[N, 0:N+1]
        A[2*(N+1)-1, 2*(N+1)-1] = A[2*(N+1)-1, 2*(N+1)-1] + k
        B[2*(N+1)-1, :] = 0.0

        # BC: real[dw_dy(xi=-1) = -k u(xi=-1)] + imag[dw_dy(xi=-1) = -k u(xi=-1)]       free stream
        A[3*(N+1)-1, :] = 0.0
        A[3*(N+1)-1, 2*(N+1):3*(N+1)] = D1[N, 0:N+1]
        A[3*(N+1)-1, 3*(N+1)-1] = A[3*(N+1)-1, 3*(N+1)-1] + k
        B[3*(N+1)-1, :] = 0.0

        # BC: real[d2pi_dxj2(xi=-1)] + imag [d2pi_dxj2(xi=-1)] = 0                      freestream      
        A[4*(N+1)-1, :] = 0.0
        A[4*(N+1)-1, 3*(N+1):4*(N+1) ] = D2[N, 0:N+1]
        A[4*(N+1)-1, 4*(N+1)-1] = A[4*(N+1)-1, 4*(N+1)-1] - k2
        B[4*(N+1)-1, :] = 0.0

        return

#============================================================================
def apply_vorticity_boundary_conditions_to_bottom_freestream(kx,kz,D1,D2,A,B):
	print 'Code not written in function apply_vorticity_boundary_conditions_to_bottom_freestream'
	print 'Rerun code using the velocity equations\n'
	sys.exit()
	return

#============================================================================
def apply_vorticity_boundary_conditions_to_top_freestream(kx,kz,D1,D2,A,B):
	print 'code not written in function apply_vorticity_boundary_conditions_to_top_freestream'
	print 'Rerun code using the velocity equations\n'
	sys.exit()
	return

#============================================================================
