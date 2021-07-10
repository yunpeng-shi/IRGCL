% set model parameters for non-uniform corruption
n=100; p=1; d=10; p_node_crpt=0.4; p_edge_crpt=0.6; crpt_type='local-adv'; 
[AdjMat, CrptMat, Ind, X, X_orig] = Nonuniform_Topology(n,d,p,p_node_crpt,p_edge_crpt,crpt_type);

% for uniform corruption use the following code
%n=100; p=1; d=10; q=0.92; crpt_type='uniform';
%[AdjMat, CrptMat, Ind, X, X_orig] = Uniform_Topology(n,d,p,q,crpt_type);


CEMP_parameters.beta_init = 1;
CEMP_parameters.beta_max = 40;
CEMP_parameters.rate = 1.2;


% IRGCL paramters
IRGCL_options.max_iter = 100;
IRGCL_options.alpha_init = 1.2;
IRGCL_options.alpha_max = 40;
IRGCL_options.rate = 1.2;
IRGCL_options.cycle_info_ratio = 1-1./((1:IRGCL_options.max_iter)+1);
IRGCL_options.LS_solver = 'PPM';
verbose = true;
[X_est, ~] = IRGCL(AdjMat, X, CEMP_parameters, IRGCL_options, verbose);
err_full = 1-sum(X_est.*X_orig,'all')/sum(X_orig,'all');
crptMat_kron = kron(CrptMat,ones(d));
err_crpt=1-sum(X_est.*X_orig.*crptMat_kron,'all')/sum(X_orig.*crptMat_kron,'all');
fprintf('mean error on all edges is %f \n', err_full);
fprintf('mean error on corrupted edges is %f \n', err_crpt);
