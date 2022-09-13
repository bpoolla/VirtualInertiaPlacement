This repository contains the MATLAB codes used in the paper:
B. K. Poolla, S. Bolognani, and F. Dörfler, "Optimal Placement of Virtual Inertia in Power Grids", IEEE Transactions on Automatic Control.

optimalplacement.m contains the raw network data: inertia, damping, network laplacian, and the optimization routine for solving the placement problem.

grad.m contains the routine which solves the rank-deficient Lyapunov solution implicitly used for solving the placement problem.

% This source code is distributed in the hope that it will be useful, but without any warranty.
% We do request that publications in which this code or modifications thereof are adopted, explicitly acknowledge that fact by citing the above mentioned paper.
