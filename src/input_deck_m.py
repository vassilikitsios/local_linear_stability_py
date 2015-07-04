# module input_deck_m.py

import math
import string
import sys

#=======================================================================
def read_deck(filename):

	file = open(filename,'r')

	EOF = False
	while not EOF:
		linestring = file.readline()
		linelist = string.split(linestring)
		print linelist
		if len(linestring) == 0:
			EOF = True
		else:
			if (linelist[0] == "equations"):
				equations = linelist[2]
				if ( (equations!="velocity_formulation") and (equations!="vorticity_formulation") ):
					print '\nFAILURE: equations != velocity_formulation or vorticity_formulation'
					sys.exit()					
			if (linelist[0] == "N"):
				N = int(linelist[2])
				if (N<=0):
					print '\nFAILURE: N<=0'
					sys.exit()					
			if (linelist[0] == "base_flow"):
				base_flow = linelist[2]
				if ( (base_flow!="couette") and (base_flow!="poiseuille") and (base_flow!="turbulent_channel") and \
				(base_flow!="blasius") and (base_flow!="shear_layer") and (base_flow!="wake") and \
				(base_flow!="file_ascii_channel_DNS") and (base_flow!="file_hdf5_TBL_DNS") and \
				(base_flow!="file_u") and (base_flow!="file_u_y") and (base_flow!="file_u_yy") ):
					print '\nFAILURE: Invalid base_flow'
					sys.exit()					
			if (linelist[0] == "base_flow_file"):
				base_flow_file = linelist[2]
			if (linelist[0] == "tol"):
				tol = float(linelist[2])
				if (tol<=0):
					print '\nFAILURE: tol<=0'
					sys.exit()					
			if (linelist[0] == "max_iter"):
				max_iter = int(linelist[2])
				if (max_iter<=0):
					print '\nFAILURE: max_iter<=0'
					sys.exit()					
			if (linelist[0] == "min_iter"):
				min_iter = int(linelist[2])
				if (min_iter<=0):
					print '\nFAILURE: min_iter<=0'
					sys.exit()					
			if (linelist[0] == "domain_size"):
				domain_size = float(linelist[2])
				if (domain_size<=0):
					print '\nFAILURE: domain_size<=0'
					sys.exit()					
			if (linelist[0] == "shear_layer_conv_u"):
				shear_layer_conv_u = float(linelist[2])
			if (linelist[0] == "u_deficit"):
				u_deficit = float(linelist[2])
			if (linelist[0] == "wake_half_width"):
				wake_half_width = float(linelist[2])
				if (wake_half_width<=0):
					print '\nFAILURE: wake_half_width<=0'
					sys.exit()					
			if (linelist[0] == "top_boundary_condition"):
				top_boundary_condition = linelist[2]
				if ( (top_boundary_condition!="wall") and (top_boundary_condition!="freestream") ):
					print '\nFAILURE: top_boundary_condition != wall or freestream'
					sys.exit()					
			if (linelist[0] == "bottom_boundary_condition"):
				bottom_boundary_condition = linelist[2]
				if ( (bottom_boundary_condition!="wall") and (bottom_boundary_condition!="freestream") ):
					print '\nFAILURE: bottom_boundary_condition != wall or freestream'
					sys.exit()					
			if (linelist[0] == "Re"):
				Re = float(linelist[2])
				if (Re<=0):
					print '\nFAILURE: Re<=0'
					sys.exit()					
                        if (linelist[0] == "kx_real_min"):
				kx_real_min = float(linelist[2])
                        if (linelist[0] == "kx_real_max"):
				kx_real_max = float(linelist[2])
                        if (linelist[0] == "N_kx_real"):
				N_kx_real = int(linelist[2])
				if (N_kx_real<=0):
					print '\nFAILURE: N_kx_real<=0'
					sys.exit()					
                        if (linelist[0] == "kx_imag_min"):
				kx_imag_min = float(linelist[2])
                        if (linelist[0] == "kx_imag_max"):
				kx_imag_max = float(linelist[2])
                        if (linelist[0] == "N_kx_imag"):
				N_kx_imag = int(linelist[2])
				if (N_kx_imag<=0):
					print '\nFAILURE: N_kx_imag<=0'
					sys.exit()					
                        if (linelist[0] == "N_kz_real"):
				N_kz_real = int(linelist[2])
				if (N_kz_real<=0):
					print '\nFAILURE: N_kz_real<=0'
					sys.exit()					
                        if (linelist[0] == "kz_real_max"):
				kz_real_max = float(linelist[2])
                        if (linelist[0] == "kz_real_min"):
				kz_real_min = float(linelist[2])
                        if (linelist[0] == "kz_imag"):
				kz_imag = float(linelist[2])
                        if (linelist[0] == "scale_domain"):
				scale_domain = linelist[2]
                        if (linelist[0] == "half_grid_point"):
				half_grid_point = float(linelist[2])
                        if (linelist[0] == "num_modes_to_write"):
				num_modes_to_write = int(linelist[2])
                        if (linelist[0] == "run_type"):
				run_type = linelist[2]
				if ( (run_type!="write_grid") and (run_type!="stability_analysis") \
					and (run_type!="group_modes") and (run_type!="transient_growth") ):
					print '\nFAILURE: Invalid run_type'
					sys.exit()					
                        if (linelist[0] == "dt"):
				dt = float(linelist[2])
                        if (linelist[0] == "num_time_steps"):
				num_time_steps = int(linelist[2])
                        if (linelist[0] == "x_pos0"):
				x_pos0 = int(linelist[2])
                        if (linelist[0] == "x_posN"):
				x_posN = int(linelist[2])

	file.close()

	return run_type, equations, N, scale_domain, base_flow, base_flow_file, tol, max_iter, min_iter,\
		domain_size, half_grid_point, shear_layer_conv_u, u_deficit, wake_half_width,\
		top_boundary_condition, bottom_boundary_condition,\
		Re, kx_real_min, kx_real_max, N_kx_real, kx_imag_min, kx_imag_max, N_kx_imag,\
		N_kz_real, kz_real_max, kz_real_min, kz_imag, num_modes_to_write, dt, num_time_steps, x_pos0, x_posN

