%%%%%%%% Optimal Placement of Inertia %%%%%%%%%

%Our paper on Inertia Placement: B. K. Poolla, S. Bolognani, and F. Dörfler, "Optimal Placement 
%of Virtual Inertia in Power Grids", IEEE Transactions on Automatic Control.

m0=[37.2423; 12.4141; 7.2193; 35.3801; 35.3801; 12.7324; 37.2423; 37.2423; 19.9899]; %vector of existing inertia allocations
Mmin= m0;
Mmax= 4*m0; %vector of maximum allowable inertia allocations
Mbdg= 405; %total inertia budget

d0= [1.9465; 3.3423; 2.7072; 2.2886; 1.1141; 4.7746; 2.2282; 2.2282; 7.4962]; %vector of existing damping

L=[ 20.2184   -7.8189  -12.3995        0         0         0         0         0         0
   -7.8189   38.8168  -30.9979         0         0         0         0         0         0
  -12.3995  -30.9979   58.8944         0         0   -7.8123         0         0   -7.6848
         0         0         0   20.5508   -7.8096  -12.7412         0         0         0
         0         0         0   -7.8096   40.4839  -32.6743         0         0         0
         0         0   -7.8123  -12.7412  -32.6743   61.6735         0         0   -8.4457
         0         0         0         0         0         0   19.7916   -7.8049  -11.9867
         0         0         0         0         0         0   -7.8049   38.3717  -30.5668
         0         0   -7.6848         0         0   -8.4457  -11.9867  -30.5668   58.6840]; %network Laplacian

options = optimoptions('fmincon', 'Display', 'iter', 'GradObj', 'off', 'GradCons', 'off', 'MaxFunctionEvaluations', 10000);

%For \ell-1 norm minimization: m_0 is initial inertia, k is the weighting factor
%Use the function [y,deltay,P0,P1]= DualGrad(m,d,L,k,m_0)

%Relative weights on penalty: w is the weight
%Use the function [y,deltay, P0, P1]= DualGrad(m,d,L,w)

%Individual capacity constraint
[xopt,fval,exitflag,output,lambda,grad]= fmincon(@(x) grad(x,d0,-L), m0, [], [], [], [], m0, Mmax, [], options);

%Total budget constant
%[xopt,fval,exitflag,output,lambda,grad]= fmincon(@(x) DualGrad(x,d0,-L), m0, ones(1,9), Mbdg, [], [], m0, [], [], options);

H2opt=fval;
Mopt=xopt;