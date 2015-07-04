# module baseflow_m.py

import math
import numpy as np
import string
import tables

import io_m as io

#============================================================================
def generate_base_flow(base_flow,scale_domain,tol,max_iter,min_iter,shear_layer_conv_u,wake_half_width,u_deficit,base_flow_file,\
			domain_size,half_grid_point,x_pos0,x_posN,eta,D1,D2,D3,D4,u,u_y,u_yy,Re):
        if (base_flow == "couette"):
                generate_couette_profile(eta,u,u_y,u_yy)
        elif (base_flow == "poiseuille"):
                generate_poiseuille_profile(eta,u,u_y,u_yy)
        elif (base_flow == "turbulent_channel"):
                generate_turbulent_channel_profile(eta,u,u_y,u_yy)
        elif (base_flow == "blasius"):
                generate_blasius_profile(tol,max_iter,min_iter,eta,D1,D2,D3,u,u_y,u_yy)
        elif (base_flow == "shear_layer"):
                generate_shear_layer_profile(shear_layer_conv_u,eta,u,u_y,u_yy)
        elif (base_flow == "wake"):
                generate_wake_profile(wake_half_width,u_deficit,eta,u,u_y,u_yy)
        elif (base_flow == "file_u"):
                read_from_file(1,base_flow_file,eta,u,u_y,u_yy)
        elif (base_flow == "file_u_y"):
               	read_from_file(2,base_flow_file,eta,u,u_y,u_yy)
        elif (base_flow == "file_u_yy"):
               	read_from_file(3,base_flow_file,eta,u,u_y,u_yy)
        elif (base_flow == "file_ascii_channel_DNS"):
                read_ascii_channel_DNS_from_file(base_flow_file,eta,u,u_y,u_yy)
        elif (base_flow == "file_hdf5_TBL_DNS"):
                read_hdf5_TBL_DNS_from_file(base_flow_file,x_pos0,x_posN,eta,u,u_y,u_yy,Re)
	io.write_base_flow(eta,u,u_y,u_yy)
	io.plot_base_flow(eta,u,u_y,u_yy)
	return

#============================================================================
def read_hdf5_TBL_DNS_from_file(filename,x_pos0,x_posN,eta,u,u_y,u_yy,Re):
        print '\n   Reading ', filename, ' ...'
	input_file = tables.openFile(filename, 'r')
	y_read = input_file.root.y[1:-1]
	N_read = len(y_read)-1
	N = len(eta)-1

	u_read = np.array(np.zeros((N_read+1), dtype=np.float))
	u0 = np.array(np.zeros((N+1), dtype=np.float))
	u_y0 = np.array(np.zeros((N+1), dtype=np.float))
	u_yy0 = np.array(np.zeros((N+1), dtype=np.float))

	num_profiles = 0
	delta99_avg = 0.0
	u_inf_avg = 0.0
        for i in range(x_pos0,x_posN+1):
		print "      processing profile ", i, " from ", x_pos0, " to ", x_posN
		u_read = input_file.root.ua[i,:]
		u_inf = np.max(u_read)
		delta99 = y_read[np.argmax(u_read)]
		print "         delta99 = ", delta99
		if (eta[-1]<y_read[-1]/delta99):
			#print u_inf, np.argmax(u_read), u_read[np.argmax(u_read)], delta99
			interpolate_onto_new_grid_using_cubic_spline(N_read+1,y_read/delta99,u_read/u_inf,eta,u0,u_y0,u_yy0)
			u += u0
			u_y += u_y0
			u_yy += u_yy0
			delta99_avg += delta99
			u_inf_avg += u_inf
			num_profiles += 1
		else:
			print "         stability grid outside of flow field grid, not included in interpolation."
	u /= num_profiles
	u_y /= num_profiles
	u_yy /= num_profiles
	delta99_avg /= num_profiles
	u_inf_avg /= num_profiles
	#interpolate_onto_new_grid_using_cubic_spline(len(eta),eta,u,eta,u,u_y,u_yy)

	input_file.close()

	u_inf = np.max(u)
	delta99_pos = np.argmax(u)
	delta99 = eta[delta99_pos]
	delta = np.sum((1.0-u[0:delta99_pos]/u_inf)*(eta[1:delta99_pos+1]-eta[0:delta99_pos]))
	theta = np.sum((1.0-u[0:delta99_pos]/u_inf)*u[0:delta99_pos]/u_inf*(eta[1:delta99_pos+1]-eta[0:delta99_pos]))

	print "\n   Averaged over ", num_profiles, " profiles:"
	print "      freestream velocity = ", u_inf_avg
	print "      boundary layer thickness = ", delta99_avg

	print "\n   Stability scaled:"
	print "      freestream velocity from profiles = ", u_inf_avg
	print "      boundary layer thickness = ", delta99
	print "      displacement thickness = ", delta
	print "      momentum thickness = ", theta
	print "      shape factor = ", delta/theta

	print "\n   Original Reynolds number = ", Re
	Re = Re*u_inf_avg*delta99_avg
	print "\n   Stability scaled Reynolds number = ", Re
	return

#============================================================================
def read_ascii_channel_DNS_from_file(filename,eta,u,u_y,u_yy):
	N = len(eta)-1
	N_read_temp = 10*N
	y_read = np.array(np.zeros((N_read_temp), dtype=np.float))
	u_read = np.array(np.zeros((N_read_temp), dtype=np.float))

        print '\n   Reading ', filename, ' ...'
        file = open(filename,'r')
        for i in range(0,27):
                linestring = file.readline()
	i = 0
	while 1:
                linestring = file.readline()
		if not linestring: break
                linelist = string.split(linestring)
                y_read[i]= float(linelist[0])-1.0
                u_read[i]= float(linelist[2])
		i = i + 1
        file.close()
	N_read = 2*(i-1)

        for i in range(N_read/2,N_read+1):
		y_read[i] = -y_read[N_read-i]
		u_read[i] = u_read[N_read-i]
	u_read[:] = u_read[:] / u_read[N_read/2]
	print "Non-dimensionalised velocity by freestrem value of: ", u_read[N_read/2]

	if (N_read!=N):
		#interpolate_onto_new_grid(N_read,y_read,u_read,eta,u)
		interpolate_onto_new_grid_using_cubic_spline(N_read+1,y_read,u_read,eta,u,u_y,u_yy)
	else:
		u_y[:] = np.dot(u,D1)
		u_yy[:] = np.dot(u,D2)
	return

#============================================================================
def read_from_file(num_fields,filename,eta,u,u_y,u_yy):
	N = len(eta)-1
	N_read_temp = 10*N
	y_read = np.array(np.zeros((N_read_temp), dtype=np.float))
	u_read = np.array(np.zeros((N_read_temp), dtype=np.float))

        print '\n   Reading ', filename, ' ...'
        file = open(filename,'r')
        for i in range(0,3):
                linestring = file.readline()
	i = 0
	while 1:
                linestring = file.readline()
		if not linestring: break
                linelist = string.split(linestring)
                y_read[i]= float(linelist[0])
                u_read[i]= float(linelist[1])
		if (num_fields>1):
                	u_y_read[i]= float(linelist[2])
		if (num_fields>2):
                	u_yy_read[i]= float(linelist[3])
		i = i + 1
	N_read = i-1
        file.close()

	if (N_read!=N):
		interpolate_onto_new_grid(N_read,y_read,u_read,eta,u)
		if (num_fields>1):
			interpolate_onto_new_grid(N_read,y_read,u_y_read,eta,u_y)
		if (num_fields>2):
			interpolate_onto_new_grid(N_read,y_read,u_yy_read,eta,u_yy)
	if (num_fields<2):
		u_y[:] = np.dot(u,D1)
	if (num_fields<3):
		u_yy[:] = np.dot(u_y,D1)
	return

#============================================================================
def interpolate_onto_new_grid(N_read,y_read,u_read,eta,u):
	N = len(eta) - 1
        u[0] = u_read[0]
	j = 0
	print "y_read_max=",y_read[-1]
	print "y_max=",eta[-1]
        for i in range(1,N+1):
		while (eta[i]>y_read[j]):
			j = j + 1
                m = (u_read[j] - u_read[j-1]) / (y_read[j] - y_read[j-1])
        	u[i] = u_read[j-1] + m*(eta[i] - y_read[j-1])
        #u[N] = u_read[N_read]
	return

#============================================================================
def interpolate_onto_new_grid_using_cubic_spline(N_read,y_read,u_read,eta,u,u_y,u_yy):
	N = len(eta)
        A = np.array(np.zeros((N_read), dtype=np.float64))
        B = np.array(np.zeros((N_read-1), dtype=np.float64))
        C = np.array(np.zeros((N_read), dtype=np.float64))
        D = np.array(np.zeros((N_read-1), dtype=np.float64))

        Q = np.array(np.zeros((N_read), dtype=np.float64))
        dy = np.array(np.zeros((N_read), dtype=np.float64))
        diag1 = np.array(np.zeros((N_read-1), dtype=np.float64))
        diag2 = np.array(np.zeros((N_read), dtype=np.float64))
        diag3 = np.array(np.zeros((N_read-1), dtype=np.float64))

        dy[0:N_read-1] = y_read[1:N_read]-y_read[0:N_read-1]
	dy[N_read-1] = dy[N_read-2]

        # assign A coefficients
        A[:] = u_read[:]

        diag1[:] = dy[0:N_read-1]
        diag1[N_read-2] = 0.0
        diag2[:] = 4.0*dy[:]
        diag2[0] = 1.0
        diag2[N_read-1] = 1.0
        diag3[:] = dy[1:N_read]
        diag3[0] = 0.0
        Q[0] = 0.0
        for i in range(1,N_read-1):
        	Q[i] = 3.0/dy[i]*(A[i+1] - 2.0*A[i] + A[i-1])
        Q[N_read-1] = 0.0

        for i in range(1,N_read-1):
                diag2[i] = diag2[i] - diag3[i-1] * diag1[i-1] / diag2[i-1]
                Q[i] = Q[i] - Q[i-1] * diag1[i-1] / diag2[i-1]

        C[N_read-1] = 0.0
        for i in range(N_read-2,0,-1):
                C[i] = (Q[i] - diag3[i] * C[i+1]) / diag2[i]

        for i in range(0,N_read-1):
                D[i] = (C[i+1] - C[i]) / 3.0 / dy[i]

        for i in range(0,N_read-1):
                B[i] = (A[i+1] - A[i]) / dy[i] - dy[i]/3.0*(2.0*C[i] + C[i+1])

	#print eta
	#print y_read
	#print len(eta)
	#print len(y_read)
        u[0] = u_read[0]
	j = 0
        for i in range(1,N):
		while (eta[i]>y_read[j]):
			j = j + 1
		j = j -1
                dy_i    = eta[i] - y_read[j]
        	u[i]    = A[j] + B[j]*dy_i +     C[j]*math.pow(dy_i,2.0) +     D[j]*math.pow(dy_i,3.0)
        	u_y[i]  =        B[j]      + 2.0*C[j]*dy_i               + 3.0*D[j]*math.pow(dy_i,2.0)
        	u_yy[i] =                  + 2.0*C[j]                    + 6.0*D[j]*dy_i
        #u[N-1] = u_read[N_read-1]
	return

#============================================================================
def generate_wake_profile(wake_half_width,u_deficit,eta,u,u_y,u_yy):
	N = len(u)-1
	for i in range(0,N+1):
		sech_n = 1.0 / np.cosh(eta[ii] / wake_half_width)
		u[i] = 1.0 - u_deficit * pow(sech_n,2.0)
		u_y[ii] = 2.0 * u_deficit * np.pow(sech_n,2.0) * np.tanh(eta[i])
		u_yy[ii] = 0.0
	print "############# Need to calculate analytical function for u_yy ###############"
	return 

#============================================================================
def generate_couette_profile(eta,u,u_y,u_yy):
	N = len(u)-1
	for i in range(0,N+1):
		u[i] = eta[i]
		u_y[i] = 1.0
		u_yy[i] = 0.0 
	return

#============================================================================
def generate_turbulent_channel_profile(eta,u,u_y,u_yy):
	N = len(u)-1
	for i in range(0,N+1):
		u[i] = 1.14*(2.01/np.pi/4.0*np.log(np.abs(np.cos(eta[ii]*np.pi/2.01))) / 7.7603807704e-01 + 1.0)
		u_y[i] = -np.tan(eta[i]*np.pi/2.01)/4.0 / 7.7603807704e-01 * 1.14
		u_yy[i] = 0.0
	u[0] = 0.0
	u[N] = 0.0
	print "############# Need to calculate analytical function for u_yy ###############"
	return

#============================================================================
def generate_poiseuille_profile(eta,u,u_y,u_yy):
	N = len(u)-1
	for i in range(0,N+1):
		u[i] = 1.0 - pow(eta[i],2.0)
		u_y[i] = -2.0 * eta[i]
		u_yy[i] = -2.0 
	return

#============================================================================
def generate_shear_layer_profile(conv_u,eta,u,u_y,u_yy):
	N = len(u)-1
	for i in range(0,N+1):
		u[i] = 1.0 + conv_u * np.tanh(eta[i])
		u_y[i] = conv_u * ( 1.0 - np.pow(np.tanh(eta[i]), 2.0) )
		u_yy[i] = 0.0
	print "############# Need to calculate analytical function for u_yy ###############"
	return

#============================================================================
def generate_blasius_profile(tol,max_iter,min_iter,eta,D1,D2,D3,u,u_y,u_yy):
	N = len(u)-1
	f = np.array(np.zeros((N+1), dtype=np.float))
	for i in range(0,N+1):
		f[i] = eta[i] - 1.21 * (1.0 - np.exp(-1.0*eta[i]))
	u[:] = np.dot(D1,f)
	u_y[:] = np.dot(D2,f)
	u_yy[:] = np.dot(D3,f)
	print "############# Note approximate profile, not iterated to get true Blasius ###############"
	return

#============================================================================
