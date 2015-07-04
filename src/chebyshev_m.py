# module chebyshev_m.py

import math
import numpy as np

#============================================================================
def scale_domain_and_derivatives(scale_domain,domain_size,half_grid_point,eta,D1,D2,D3,D4):
        if (scale_domain == "finite"):
                scale_for_finite_domain(domain_size,eta,D1,D2,D3,D4)
        elif(scale_domain == "finite_unclustered"):
                scale_for_unclustered_finite_domain(domain_size,eta,D1,D2,D3,D4)
        elif(scale_domain == "semi_finite"):
                scale_for_semi_finite_domain(domain_size,eta,D1,D2,D3,D4)
        elif(scale_domain == "semi_finite_unclustered"):
                scale_for_unclustered_semi_finite_domain(domain_size,eta,D1,D2,D3,D4)
        elif(scale_domain == "semi_finite_clustered"):
                scale_for_clustered_semi_finite_domain(domain_size,half_grid_point,eta,D1,D2,D3,D4)

#============================================================================
def generate_chebyshev_integral_weights(N,IW):
        if (N%2==0):	 # EVEN N
                IW[0] = 1.0 / (1.0 * math.pow(1.0*N,2.0) - 1.0)
        	for i in range(1,N):
        		for k in range(0,N/2+1):
                                c = 1.0
                                if (k==0):
					c = 2.0
                                IW[i] = IW[i] + 4.0 / N / c * np.cos(2.0*math.pi*i*k/N) / (1.0-4.0*math.pow(1.0*k,2.0))
                IW[N] = IW[0]
        else:		 # ODD N
                IW[0] = 1.0 / (1.0 * math.pow(1.0*N,2.0))
        	for i in range(1,N):
        		for k in range(0,(N-1)/2+1):
                                c = 1.0
                                if ( (k==0) or (k==(N-1)/2) ):
					c = 2.0
                                IW[i] = IW[i] + 4.0 / N / c * np.cos(2.0*math.pi*i*k/N) / (1.0-4.0*math.pow(1.0*k,2.0)) \
                                	+ 4.0 / N * pow(-1.0,i) * 2.0 * np.cos(math.pi*i/N) / (1.0*math.pow(1.0*N,2.0)*(2.0-N*1.0))
                IW[N] = IW[0]
	return

#======================================================================-
def generate_grid_and_derivative_matrices(eta,D1,D2,D3,D4):
	generate_colocation_points(eta)
	generate_D1(D1,eta)
	generate_D2(D2,eta)
	D3[:,:] = np.dot(D1,D2)
	correct_diagonal_terms(D3)
	D4[:,:] = np.dot(D2,D2)
	correct_diagonal_terms(D4)
	#test_all_derivative_matrices(eta,D1,D2,D3,D4)
	return 

#======================================================================-
def scale_for_finite_domain(domain_size,eta,D1,D2,D3,D4):
	eta[:] = -eta[:] / 2.0 * domain_size
	scale_derivatives_for_finite_domain(domain_size,D1,D2,D3,D4)
	return 

#======================================================================-
def scale_for_semi_finite_domain(domain_size,eta,D1,D2,D3,D4):
	N0=len(eta)-1
	for i in range(0,N0+1):
		eta[i] = (1.0-eta[i]) / 2.0 * domain_size
	scale_derivatives_for_finite_domain(domain_size,D1,D2,D3,D4)
	return 

#======================================================================-
def scale_derivatives_for_finite_domain(domain_size,D1,D2,D3,D4):
	D1[:,:] = -2.0 / domain_size * D1[:,:]
        D2[:,:] = 4.0 / pow(domain_size,2.0) * D2[:,:]
        D3[:,:] = -8.0 / pow(domain_size,3.0) * D3[:,:]
        D4[:,:] = 16.0 / pow(domain_size,4.0) * D4[:,:]
	return 

#============================================================================
def scale_for_unclustered_finite_domain(domain_size,eta,D1,D2,D3,D4):
	N0=len(eta)-1
        a = cos(8*math.pi/N0)
        for i in range(0,N0+1):
                eta[i] = -domain_size * asin(a * xi[i]) / asin(a) / 2.0
        scale_derivatives_for_unclustered_finite_domain(domain_size,D1,D2,D3,D4)
        return

#============================================================================
def scale_for_shifted_unclustered_finite_domain(y_min,y_max,eta,D1,D2,D3,D4):
	N0=len(eta)-1
        a = cos(8*math.pi/N0)
        for i in range(0,N0+1):
                eta[i] = (y_min - y_max) * asin(a * xi[i]) / asin(a) / 2.0  + (y_max + y_min)/2.0
        scale_derivatives_for_unclustered_finite_domain(domain_size,eta,D1,D2,D3,D4)
        return

#============================================================================
def scale_derivatives_for_unclustered_finite_domain(domain_size,eta,D1,D2,D3,D4):
	N0=len(eta)-1
        a = cos(8*math.pi/N0)
        b = asin(a)
        for i in range(0,N0+1):
                d_xi_d_eta = -2.0 * b / domain_size / a * cos(-b * 2.0 * eta[i] / domain_size)
                d2_xi_d_eta2 = -pow(2.0 * b / domain_size,2.0) / a * sin(-2.0 * b * eta[i] / domain_size)
                for j in range(0,N0+1):
                        D2[i,j] = D2[i,j] * pow(d_xi_d_eta,2.0) + D1[i,j] * d2_xi_d_eta2
        correct_diagonal_terms(D2)

        for i in range(0,N0+1):
                d_xi_d_eta = -2.0 * b / domain_size / a * cos(b*(-2.0 * eta[i]/domain_size))
                for j in range(0,N0+1):
                        D1[i,j] = D1[i,j] * d_xi_d_eta
        correct_diagonal_terms(D1)

        D3[:,:] = np.dot(D1,D2)
        correct_diagonal_terms(D3)
        D4[:,:] = np.dot(D2,D2)
        correct_diagonal_terms(D4)
        return

#============================================================================
def scale_for_unclustered_semi_finite_domain(domain_size,eta,D1,D2,D3,D4):
	N0 = len(eta)-1
        a = math.cos(8*math.pi/N0)
        for i in range(0,N0+1):
                eta[i] = domain_size * (1.0 - math.asin(a * eta[i]) / math.asin(a)) / 2.0
        scale_derivatives_for_unclustered_semi_finite_domain(domain_size,eta,D1,D2,D3,D4)
        return

#============================================================================
def scale_derivatives_for_unclustered_semi_finite_domain(domain_size,eta,D1,D2,D3,D4):
        N0 = len(eta)-1
	a = math.cos(8*math.pi/N0)
        b = math.asin(a)
        for i in range(0,N0+1):
                d_xi_d_eta = -2.0 * b / domain_size / a * math.cos(b*(1.0 - 2.0 * eta[i]/domain_size))
                d2_xi_d_eta2 = -pow(2.0 * b / domain_size,2.0) / a * math.sin(b*(1.0 - 2.0 * eta[i]/domain_size))
                for j in range(0,N0+1):
                        D2[i,j] = D2[i,j] * pow(d_xi_d_eta,2.0) + D1[i,j] * d2_xi_d_eta2
        correct_diagonal_terms(D2)

        for i in range(0,N0+1):
                d_xi_d_eta = -2.0 * b / domain_size / a * math.cos(b*(1.0 - 2.0 * eta[i]/domain_size))
                for j in range(0,N0+1):
                        D1[i,j] = D1[i,j] * d_xi_d_eta
        correct_diagonal_terms(D1)

        D3[:,:] = np.dot(D1,D2)
        correct_diagonal_terms(D3)
        D4[:,:] = np.dot(D2,D2)
        correct_diagonal_terms(D4)
        return

#============================================================================
def scale_for_clustered_semi_finite_domain(domain_size, half_grid_point, eta, D1, D2, D3, D4):
        N0 = len(eta)-1
        denom = 1.0 - 2.0 * half_grid_point / domain_size
        if (abs(denom)<1.0e-6):
                print "cannot set half_grid_point at the mid point of the domain"
                sys.exit()
        l = half_grid_point / denom
        s = 2.0 * l / domain_size
        print "Chebyshev scaling of clustered semi infinite domain: s = ", s, ", l = ", l

        for i in range(0,N0+1):
                eta[i] = l * (1.0 - eta[i]) / (1.0 + s + eta[i])

        scale_derivatives_for_clustered_semi_finite_domain(s,l,eta,D1,D2,D3,D4)
        return

#============================================================================
def scale_derivatives_for_clustered_semi_finite_domain(s,l,eta,D1,D2,D3,D4):
        N0 = len(eta)-1
        for i in range(0,N0+1):
                d_xi_d_eta = -pow(1.0 + s + xi[i],2.0) / l / (2.0 + s)
                d2_xi_d_eta2 = 2.0 * pow(1.0 + s + xi[i],3.0) / pow(l,2.0) / pow(2.0 + s,2.0)
                for j in range(0,N0+1):
                        D2[i,j] = D2[i,j] * pow(d_xi_d_eta,2.0) + D1[i,j] * d2_xi_d_eta2
        correct_diagonal_terms(D2)

        for i in range(0,N0+1):
                d_xi_d_eta = -pow(1.0 + s + xi[i],2.0) / l / (2.0 + s)
                for j in range(0,N0+1):
                        D1[i,j] = D1[i,j] * d_xi_d_eta
        correct_diagonal_terms(D1)

        D3[:,:] = np.dot(D1,D2)
        correct_diagonal_terms(D3)
        D4[:,:] = np.dot(D2,D2)
        correct_diagonal_terms(D4)
        return

#======================================================================-
def test_derivative_matrices(eta,D1,D2):
	N0 = len(eta)-1
	test = np.array(np.zeros((N0+1,1), dtype=np.float))
	for i in range(0,N0+1):
		test[i,0] = math.pow(eta[i],2.0)
	print test
	print np.dot(D1,test)
	print np.dot(D2,test)
	return

#======================================================================-
def test_all_derivative_matrices(eta,D1,D2,D3,D4):
	N0 = len(eta)-1
	test = np.array(np.zeros((N0+1,1), dtype=np.float))
	for i in range(0,N0+1):
		test[i,0] = math.pow(eta[i],4.0)
	print test
	print np.dot(D1,test)
	print np.dot(D2,test)
	print np.dot(D3,test)
	print np.dot(D4,test)
	return

#======================================================================-
def generate_colocation_points(eta):
	N0 = len(eta)-1
	eta[:] = 0.0
	for i in range(0,N0+1):
		eta[i] = np.cos(i * np.pi / N0)
	return 

#============================================================================
def generate_D1(D1,eta):
	N0 = len(eta)-1
	D1[:,:] = 0.0
	for i in range(0,N0+1):
		for j in range(0,N0+1):
			if (i!=j):
				D1[i,j]=math.pow(-1.0,1.0*(i+j))/(eta[i]-eta[j])
	D1[0,0:N0+1]  = 2.0*D1[0,0:N0+1]
	D1[N0,0:N0+1] = 2.0*D1[N0,0:N0+1]
	D1[0:N0+1,0]  = 0.5*D1[0:N0+1,0]
	D1[0:N0+1,N0] = 0.5*D1[0:N0+1,N0]
	correct_diagonal_terms(D1)
	return 

#============================================================================
def generate_D2(D2,eta): 
	N0 = len(eta)-1
	D2[:,:] = 0.0
	for i in range(1,N0):
		for j in range(0,N0+1):
			if (i!=j):
				D2[i,j]=pow(-1.0,1.0*(i+j))*(pow(eta[i],2.0)+eta[i]*eta[j]-2.0) / ( (1.0-pow(eta[i],2.0))*pow((eta[i]-eta[j]),2.0))
	D2[0:N0+1,0]=D2[0:N0+1,0]/2.0
	D2[0:N0+1,N0]=D2[0:N0+1,N0]/2.0
	
	for j in range(1,N0+1):
		D2[0,j] = 2.0/3.0*pow(-1.0,1.0*j)*((2.0*N0*N0+1.0)*(1.0-eta[j])-6.0)/pow(1.0-eta[j],2.0)
	D2[0,N0]=D2[0,N0]/2.0
	
	for j in range(0,N0):
		D2[N0,j] = 2.0/3.0*pow(-1.0,1.0*(j+N0))*((2.0*N0*N0+1.0)*(1.0+eta[j])-6)/pow(1.0+eta[j],2.0)
	D2[N0,0]=D2[N0,0]/2.0
	
	correct_diagonal_terms(D2)
	return 

#============================================================================
def correct_diagonal_terms(D): 
	N0 = len(D)-1
	for i in range(0,N0+1):
		D[i,i] = 0.0;
		D[i,i] = -np.sum(D[i,:])
	return 

#======================================================================-
