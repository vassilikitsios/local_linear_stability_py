----------------------------------------------------------------------
1) Code overview:
----------------------------------------------------------------------

Author: Vassili Kitsios

This python code performs local stability calculations for a series of analytical profiles and profiles read in from file.

This README file contains installation instructions and an overview of the code. If you are to use this code in your own projects please cite the following documents:

Kitsios, V., Cordier, L., Bonnet, J.-P., Ooi, A. & Soria, J., 2010, Development of a nonlinear eddy viscosity closure for the stability analysis of a turbulent channel flow, Journal of Fluid Mechanics, Vol. 664, pp 74-107.
http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=7928161

Kitsios, V., 2010, Recovery of fluid mechanical modes in unsteady separated flows, PhD Thesis, The University of Melbourne
https://minerva-access.unimelb.edu.au/handle/11343/35705


----------------------------------------------------------------------
2) Installation instructions:
----------------------------------------------------------------------

This code requires only the following python libraries:

	numpy

	scipy

	matplotlib

----------------------------------------------------------------------
3) List of files and directories with brief explanations: 
----------------------------------------------------------------------

src/				- directory containing python code

tests/				- directory containing Blasius boundary layer and turbulent channel test cases

tests/blasius/			- stability analysis of a laminar Blasius boundary layer, see PhD thesis appendix B

tests/channel/			- stability analysis of a channel flow, see PhD thesis chapter 6


----------------------------------------------------------------------
3) List of files and directories within the results directories: 
----------------------------------------------------------------------

run						- batch file to execute the code

local_linear_stability.in			- input deck containing parameters specifying the EVP

results/					- directory containing the eigvalue and eigenvector files

images/						- directory containing the images created from the results


----------------------------------------------------------------------
