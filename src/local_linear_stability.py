#!/usr/bin/env python

#=======================================================================
# import libraries
import numpy as np
import string
import sys
import os

import input_deck_m as input_deck 
import chebyshev_m as cheb 
import baseflow_m as bf 
import navier_stokes_m as ns 
import group_modes_m as gm 
import io_m as io

#=======================================================================
print '\n======================================================================='
print 'Running ', sys.argv[0] 
print '=======================================================================\n'

#=======================================================================
print '\nReading input deck ...'
run_type, equations, N, scale_domain, base_flow, base_flow_file, tol, max_iter, min_iter, \
	domain_size, half_grid_point, shear_layer_conv_u, u_deficit, wake_half_width,\
	top_BC, bottom_BC, Re, kx_real_min, kx_real_max, N_kx_real, kx_imag_min, kx_imag_max, N_kx_imag,\
	N_kz_real, kz_real_max, kz_real_min, kz_imag, num_modes_to_write, \
	dt, num_time_steps, x_pos0, x_posN = input_deck.read_deck("local_linear_stability.in")

if not os.path.exists('./results'):
        os.mkdir('./results')
if not os.path.exists('./images'):
        os.mkdir('./images')

#=======================================================================
print '\nAllocating memory ...'
kx = np.array(np.zeros((N_kx_real,N_kx_imag), dtype=np.complex))	# streamwise wavenumber set
kz = np.array(np.zeros((N_kz_real), dtype=np.complex))			# spanwise wavenumber set
u = np.array(np.zeros((N+1), dtype=np.float))				# baseflow velocity
u_y = np.array(np.zeros((N+1), dtype=np.float))				# baseflow velocity first wall normal derivative
u_yy = np.array(np.zeros((N+1), dtype=np.float))			# baseflow velocity second wall normal derivative
eta = np.array(np.zeros((N+1), dtype=np.float))				# chebyshev collocation grid
D1 = np.array(np.zeros((N+1,N+1), dtype=np.float))			# chebyshev 1st derivative matrix
D2 = np.array(np.zeros((N+1,N+1), dtype=np.float))			# chebyshev 2nd derivative matrix
D3 = np.array(np.zeros((N+1,N+1), dtype=np.float))			# chebyshev 3rd derivative matrix
D4 = np.array(np.zeros((N+1,N+1), dtype=np.float))			# chebyshev 4th derivative matrix
IW = np.array(np.zeros((N+1), dtype=np.float))				# chebyshev integral weights 
if (equations=="velocity_formulation"):
	Nev = 4*(N+1)							# dimension of velocity formulation eigenvalue operators
elif (equations=="vorticity_formulation"):
	Nev = 2*(N+1)							# dimension of vorticity formulation eigenvalue operators

#=======================================================================
print '\nGenerating grid and Chebyshev matrices ...'
cheb.generate_grid_and_derivative_matrices(eta,D1,D2,D3,D4)
cheb.generate_chebyshev_integral_weights(N,IW)
cheb.scale_domain_and_derivatives(scale_domain,domain_size,half_grid_point,eta,D1,D2,D3,D4)

print '\nGenerating wavenumber set ...'
ns.generate_wavenumber_sets(kx_real_min, kx_real_max, N_kx_real, kx_imag_min, kx_imag_max, N_kx_imag,\
		                  N_kz_real, kz_real_max, kz_real_min, kz_imag, kx, kz)

if (run_type=="write_grid"):
        print "\nWriting grid to file ..."
	io.write_base_flow(eta,u,u_y,u_yy)
	io.plot_base_flow(eta,u,u_y,u_yy)
        print "\nInterpolate field onto this grid and write profile in the same format."
        print "Re-run the stability code with 'run_type = stability_analysis'\n"

elif (run_type=="stability_analysis"):
	print '\nGenerating base flow ...'
	bf.generate_base_flow(base_flow,scale_domain,tol,max_iter,min_iter,shear_layer_conv_u,wake_half_width,u_deficit,base_flow_file,\
				domain_size,half_grid_point,x_pos0,x_posN,eta,D1,D2,D3,D4,u,u_y,u_yy,Re)
	print '\nBuilding linear operators and the solving eigenvalue problem ...'
	ns.run_stability_analysis(Nev,N_kx_real,N_kx_imag,N_kz_real,kx,kz,equations,\
					eta,u,u_y,u_yy,Re,D1,D2,D4,top_BC,bottom_BC,num_modes_to_write)
	if (N_kx_real*N_kx_imag*N_kz_real>1):
        	print "To group the eigenvalues into like modes re-run the stability code with 'run_type = group_modes'\n"

elif (run_type=="group_modes"):
	print '\nGrouping eigenvalues into modes ...'
	gm.group_modes(Nev,N_kx_real,N_kx_imag,N_kz_real,kx,kz,num_modes_to_write)

elif (run_type=="transient_growth"):
	print '\nGenerating base flow ...'
	bf.generate_base_flow(base_flow,scale_domain,tol,max_iter,min_iter,shear_layer_conv_u,wake_half_width,u_deficit,base_flow_file,\
				domain_size,half_grid_point,x_pos0,x_posN,eta,D1,D2,D3,D4,u,u_y,u_yy,Re)
	print '\nBuilding linear operators and the solving eigenvalue problem ...'
	ns.run_transient_growth_analysis(Nev,N_kx_real,N_kx_imag,N_kz_real,kx,kz,equations,\
					eta,u,u_y,u_yy,Re,D1,D2,D4,IW,top_BC,bottom_BC,num_modes_to_write,dt,num_time_steps)

#=======================================================================
# EOF
#=======================================================================
