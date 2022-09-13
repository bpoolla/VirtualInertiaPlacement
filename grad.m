%%%%%%% H2 norm optimization for inertia %%%%%%%%

%For \ell-1 norm minimization: m_0 is initial inertia, k is the weighting factor
%Use the function [y,deltay,P0,P1]= grad(m,d,L,k,m_0)

%Relative weights on penalty: w is the weight
%Use the function [y,deltay, P0, P1]= grad(m,d,L,w)

function [y, deltay, P0, P1] = grad(m,d,L)
    %Consider the following approximation for the trace y evaluated in m1
    %Mathematically, y(m1 + epsilon delta) = y(m1) + epsilon y1(delta) + o(epsilon^2)
    %We know how to compute y(m1) and y1(delta)

    n= length(m); % System size
    %Omega= eye(n)-ones(n)/n; % Projection matrix
    
    %Penalties on frequency
    Z= eye(n);      %For identical penalties
    %Z= diag(d);    %For primary control penalty
   
    % System Matrices
    A0= [zeros(n), eye(n); inv(diag(m))*L, -inv(diag(m))*diag(d)]; %A matrix
    B0= [zeros(n); inv(diag(m))]; %B matrix
    Q0 = [-L, zeros(n); zeros(n), Z]; %C^T C matrix for the system with weights L, Z
   
    %Disturbance weighting matrices
    
    %T0= diag(d); % penalty scales as D
    %T0= diag([0,0,1,0,0,0,0,0,0]); %Localized disturbance @ node 4 of the
    T0= diag((1/n)*ones(1,n)); %Uniform disturbance
  
    B0= B0*(T0^0.5); %Input matrix scaled with the disturbance strengths
      
    %Determining P0, solution to the Lyapunov equation
    [XX0,~,~]= svd(Q0);

    XC=XX0(:,1:end-1);
    A0tilde=XC'*A0*XC;
    Q0tilde=XC'*Q0*XC;

    %Solve for the Observability Gramian P
    P0tilde= lyap(A0tilde',Q0tilde);
    P0= XC*P0tilde*XC';
   
    P0*A0+A0'*P0+Q0; %check if this sum is zero, i.e., P0 solves the equation 

    y= trace(B0'*P0*B0);
    %y= trace(B0'*P0*B0)+k*sum(m-m_0); %if considering an \ell-1 penalty
   
    %computing the gradient i.e., eye * g, where g is the gradient
    deltay= zeros(1,n); %initialization
    del_m=diag(ones(1,n));

    for i= 1:n            
        delta = del_m(i,:); %select i-th row of eye
        
        % Perturbed matrics
        A1= [zeros(n), zeros(n); -diag(delta)*inv(diag(m))^2*L, diag(delta)*(inv(diag(m))^2)*diag(d)];
        B1= [zeros(n); -diag(delta)*inv(diag(m))^2];
        B1=B1*(T0^0.5);

        %solve Lyapunov Equation in epsilon
        Q1= A1'*P0 + P0*A1;
        Q1tilde= XC'*Q1*XC;
        P1tilde= lyap(A0tilde',Q1tilde);
        P1= XC*P1tilde*XC';
        
        P1*A0+A0'*P1+Q1; %check if this sums up to zero

        deltay(i)= trace(B1'*P0*B0 + B0'*P1*B0 + B0'*P0*B1);
        %deltay(i)= trace(B1'*P0*B0 + B0'*P1*B0 + B0'*P0*B1)+k; %for \ell-1 
    end
end


