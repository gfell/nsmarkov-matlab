function [y_grid, trans] = lcrouwenhorst(rho,sigma_eps,N,T) 
% ========================================================================
%  Code for discretization method in "An Extension of Rouwenhorst Method
%  for Non-Stationary Environments" by Giulio Fella, Giovanni Gallipoli
%  and Jutong Pan
% 
%  Rouwenhurst method to approximate non-stationary AR(1) process by a discrete Markov chain
%       y(t) = rho(t)*y(t-1)+ epsilon(t),   epsilion(t)~iid N(0,sigmaeps(t))
%       with INITIAL condition y(0) = 0 (equivalently y(1)=epsilon(1) ) 
%  INPUT:  rho - Tx1 vector of serial correlation coefficients
%          sigma_eps - Tx1 vector of standard deviations of innovations
%  OUTPUT: trans  - NxNxT array of T (NxN) transition matrices  
%                   transition probabilities are arranged by row
%          y_grid - NxT matrix, each column stores the Markov state 
% 	   space for period t
% !========================================================================%

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

% *** Step 1: construct the state space y_grid(t) in each period t.
% Evenly-spaced N-state space over [-sqrt(N-1)*sigma_y(t),sqrt(N-1)*sigma_y(t)].


% 1.a Compute unconditional variances of y(t)
sigma_y(1) = sigma_eps(1);
for i = 2:T 
    sigma_y(i) = sqrt(rho(i)^2*sigma_y(i-1)^2+sigma_eps(i)^2);
end

% 1.b Construct state space
h = 2*sqrt(N-1)*sigma_y/(N-1); % grid step
y_grid = repmat(h,N,1);
y_grid(1,:)=-sqrt(N-1) * sigma_y;
y_grid = cumsum(y_grid,1); 

%  *** Step 2: Compute the transition matrices trans(:,:,t) from
%              period (t-1) to period t
% The transition matrix for period t is defined by parameter p(t).
% p(t) = 0.5*(1+rho*sigma(t-1)/sigma(t))

% Note: trans(:,:,1) is the transition matrix from y(0)=0 to any
% gridpoint of y_grid(1) in period 1.
% Any of its rows is the (unconditional) distribution in period 1.

p = 1/2; % First period: p(1) = 0.5 as y(1) is white noise.  
trans(:,:,1) = rhmat(p,N);

for j = 2:T 
        % Compute p for t>1 
        p = (sigma_y(j)+rho(j)*sigma_y(j-1))/(2*sigma_y(j) );   
        trans(:,:,j) = rhmat(p,N);
end

    
    function [Pmat] = rhmat(p,N)
    % Computes Rouwenhorst matrix as a function of p and N
        
    Pmat = zeros(N,N);
    
    % Step 2(a): get the transition matrix P1 for the N=2 case
    if N == 2
        Pmat = [p, 1-p; 1-p, p];  
    else
        P1 = [p, 1-p; 1-p, p];
    
        % Step 2(b): if the number of states N>2, apply the Rouwenhorst
        % recursion to obtain the transition matrix trans
        for i = 2:N-1
            P2 = p *     [P1,zeros(size(P1,1),1); zeros(1,size(P1,2)),0 ] + ...
                 (1-p) * [zeros(size(P1,1),1),P1; 0,zeros(1,size(P1,2)) ] + ...
                 (1-p) * [zeros(1,size(P1,2)),0 ; P1,zeros(size(P1,1),1)] + ...
                 p *     [0,zeros(1,size(P1,2)) ; zeros(size(P1,1),1),P1];
    
            P2(2:i,:) = 0.5*P2(2:i,:);
            
            if i==N-1
                Pmat = P2;
            else
                P1 = P2;
            end
        end % of for
    end % if N == 2
    
    end % of rhmat function
end % function lcrouwenhorst