%% Author: Yunpeng Shi and Shaohan Li
%%------------------------------------------------
%% Robust Multi-object Matching via Iterative Reweighting of the Graph Connection Laplacian
%%------------------------------------------------
%% Input Parameters: 
%% AdjMat: n by n adjacency matrix
%% X: block matrix of relative permutations. Each [i,j]-th block is 
%%    the d by d relative permutation Xij between node i and j

%% CEMP_parameters.beta_init: the rweighting parameter of CEMP in the 1st iteration
%% CEMP_parameters.rate:   in each iteration, beta is increased by multiplying the rate
%% CEMP_parameters.beta_max: whenever beta exceeds beta_max, stop CEMP iteration 

%% IRGCL_options.alpha_init: the rweighting parameter of IRGCL in the 1st iteration
%% IRGCL_options.rate:   in each iteration, alpha is increased by multiplying the rate
%% IRGCL_options.alpha_max: whenever alpha exceeds alpha_max, stop IRGCL iteration 
%% IRGCL_options.max_iter: the maximal number of iterations of IRGCL
%% IRGCL_options.cycle_info_ratio: the coefficient lambda_t of cycle-consistency information.
%% IRGCL_options.LS_solver: The choice of solvers for the weighted least squares in each iteration;
%% We offer two options: 'PPM' refers to the projected power method, 'Spectral' refers to eigen-decomposition method.


%% Output:
%% X_est: nd by nd block matrix of the estimated relative permutions

%% Reference
%% [1] Yunpeng Shi, Shaohan Li and Gilad Lerman. 
%% "Robust Multi-object Matching via Iterative Reweighting of the Graph Connection Laplacian" NeurIPS 2020.



function[X_est] = IRGCL(AdjMat, X, CEMP_parameters, IRGCL_options)

% Check existence of 3-cycles
CoDeg = sign((AdjMat^2).*AdjMat);
if norm(CoDeg-AdjMat,'fro')>0
    error('Error. Every edge must be contained in at least one 3-cycle');
end

%CEMP parameters
beta_init = CEMP_parameters.beta_init;
beta_max = CEMP_parameters.beta_max;
rate_cemp = CEMP_parameters.rate;


% IRGCL paramters
maxIters = IRGCL_options.max_iter;
alpha_init = IRGCL_options.alpha_init;
alpha_max = IRGCL_options.alpha_max;
rate_irgcl = IRGCL_options.rate;
lambda_irgcl = IRGCL_options.cycle_info_ratio;
T_lambda_irgcl = length(lambda_irgcl);
if T_lambda_irgcl<maxIters
    lambda_irgcl = [lambda_irgcl,lambda_irgcl(end)*(ones(1,maxIters-T_lambda_irgcl))];
end

n = size(AdjMat,1);
N = size(X,1);
d=N/n;



Weights=AdjMat;
conv_filter = @(block_struct) mean(block_struct.data(:));


beta = beta_init;
delta = 1e-8; % try 1e-16 if the weight matrix is ill-conditioned.
while beta <= beta_max   
    % parameter controling the decay rate of reweighting function
    Weight_cemp=kron(Weights, ones(d));    
    S = X.*Weight_cemp;
    S2 = S^2;
    Weights2 = max(Weights^2, delta);
    Weight_cemp2 = kron(Weights2, ones(d));
    ProdS2X = (S2./Weight_cemp2).*X;
    A = blockproc(ProdS2X, [d d], conv_filter).*AdjMat.*d;
    Weights = exp(-beta.*(1-A)).*AdjMat;
    beta = beta*rate_cemp;
end

Weights = A.*AdjMat;
Weights = diag(1./max(sum(Weights,2), delta))*Weights;
Weight_cemp=kron(Weights, ones(d));
S = X.*Weight_cemp;
[V,~] = eigs(S,d,'la');
PMat = zeros(n*d,d);
for i=1:n
    PMat(((i*d-d+1):(i*d)), :) = project_hungarian( V(((i*d-d+1):(i*d)), :)*( V(1:d, :))');
end
X_est = PMat*PMat';

score = inf;
iter=1;
stop_threshold = 1e-10;
alpha = alpha_init;
while((score>stop_threshold)&&(iter<maxIters))
    lam = lambda_irgcl(iter);
    Prod1 = X_est.*X;
    A1 = blockproc(Prod1, [d d], conv_filter).*AdjMat.*d;
    Weights = exp(-alpha.*(1-A1)).*AdjMat;
    Weight_cemp=kron(Weights, ones(d));    
    S = X.*Weight_cemp;
    S2 = S^2;
    Weights2 = max(Weights^2,delta);
    Weight_cemp2 = kron(Weights2, ones(d));
    ProdS2X = (S2./Weight_cemp2).*X;
    A2 = blockproc(ProdS2X, [d d], conv_filter).*AdjMat.*d;
    Weights = (1-lam).*A1+lam.*A2;
   
    
    if strcmp(IRGCL_options.LS_solver,'PPM')    
        Weight_cemp=kron(Weights, ones(d));
        S = X.*Weight_cemp;
        PMat = S*PMat;
        for i=1:n
            PMat(((i*d-d+1):(i*d)), :) = project_hungarian( PMat(((i*d-d+1):(i*d)), :) );
        end
    end
    if strcmp(IRGCL_options.LS_solver,'Spectral')
        Weights = diag(1./max(sum(Weights,2),delta))*Weights;
        Weight_cemp=kron(Weights, ones(d));
        S = X.*Weight_cemp;
        [V,~] = eigs(S,d,'la');
        PMat = zeros(n*d,d);
        for i=1:n
            PMat(((i*d-d+1):(i*d)), :) = project_hungarian( V(((i*d-d+1):(i*d)), :)*( V(1:d, :))');
        end
    end
   
    
    X_est_new = PMat*PMat'; 
    score = (norm(X_est_new-X_est,'fro'))^2/(2*d);
    iter = iter+1;
    X_est = X_est_new;
    alpha = min(rate_irgcl*alpha, alpha_max);
end