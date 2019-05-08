function [y_grid, trans] = lctauchen(rho,sigma_eps,omega,N,T)    
% ========================================================================
%  Code for discretization method in "An Extension of Rouwenhorst Method
%  for Non-Stationary Environments" by Giulio Fella, Giovanni Gallipoli
%  and Jutong Pan
% 
%  Tauchen method to approximate non-stationary AR(1) process by a discrete Markov chain
%       y(t) = rho(t)*y(t-1)+ epsilon(t),   epsilion(t)~iid N(0,sigmaeps(t))
%       with INITIAL condition y(0) = 0 (equivalently y(1)=epsilon(1) ) 
%  INPUT:  rho 	     - Tx1 vector of serial correlation coefficients
%          sigma_eps - Tx1 vector of standard deviations of innovations
%          omega     - Tx1 vector of scaling factors for grid range
%                      (equal 3 in Tauchen (1986) paper)
%  OUTPUT: trans     - NxNxT matrix of T (NxN) transition matrices  
%                      transition probabilities are arranged by row
%          y_grid    - an NxT matrix, each column stores the Markov state 
% 	               space for period t
% ========================================================================%

    sigma_y = zeros(1,T);
    y_grid = zeros(N,T); 
    trans = zeros(N,N,T);

    if N < 2
        disp('The state space has to have dimension N>1. Exiting.')
        return;
    end

    if T < 2
        disp('The time horizon has to have dimension N>1. Exiting.')
        return;
    end

    % *** Step 1: construct the state space y_grid(t) in each period t,
    % Evenly-spaced N-state space over [-omega(t)*sigma_y(t),omega(t)*sigma_y(t)].

    % 1. Compute unconditional variances of y(t)
    sigma_y(1) = sigma_eps(1);
    for i = 2:T 
        sigma_y(i) = sqrt(rho(i)^2*sigma_y(i-1)^2+sigma_eps(i)^2);
    end
    
    
    % Construct state space
    h = 2*omega'.*sigma_y/(N-1); % grid step
        
    y_grid = repmat(h,N,1);
    y_grid(1,:)=-omega(1)* sigma_y;
    y_grid = cumsum(y_grid,1);

    %  *** Step 2: Compute the transition matrices trans(:,:,t) from
    %              period (t-1) to period t
    
    % Compute the transition matrix in period 1; i.e., from y(0)=0
    % to any gridpoint of y_grid(1) in period 1.
    % Any of its rows is the (unconditional) distribution in period 1.
    for i =1:N
        temp1d = (y_grid(:,1)-h(1)/2)/sigma_eps(1);
        temp1d = max(temp1d,-37); % To avoid underflow in next line
        cdf(i,:) = cdf_normal(temp1d);
    end 
    trans(:,1,1) = cdf(:,2); 
    trans(:,N,1) = 1-cdf(:,N);
    for j=2:N-1
        trans(:,j,1) = cdf(:,j+1)-cdf(:,j);
    end
 
    % Compute the transition matrices for t>2
    for t=2:T
        for i=1:N
            temp3d(i,:,t) = (y_grid(:,t) - rho(t)*y_grid(i,t-1) - h(t)/2)/sigma_eps(t);
            temp3d(i,:,t) = max(temp3d(i,:,t),-37); % To avoid underflow in next line
            cdf(i,:) = cdf_normal(temp3d(i,:,t));
        end 
        trans(:,1,t) = cdf(:,2);   
        trans(:,N,t) = 1-cdf(:,N);
        for j=2:N-1
            trans(:,j,t) = cdf(:,j+1)-cdf(:,j);
        end
    end

    function c = cdf_normal(x)
    % Returns the value of the cdf of the Standard Normal distribution at point x
        
        c = 0.5 * erfc(-x/sqrt(2));
    end % function cdf_normal
end % function lctauchen