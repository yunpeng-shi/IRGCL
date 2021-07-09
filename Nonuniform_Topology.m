%% Author: Yunpeng Shi and Shaohan Li
%%------------------------------------------------
%% function for generating data that follows non-uniform corruption model
%%------------------------------------------------
%% Input Parameters: 
%% n: the number of graph nodes (absolute permutations)
%% d: dimension of each permutation matrix (d by d)
%% p: probability of connecting two nodes (controls sparsity, p=1 refers to a complete graph)
%% p_node_crpt: probability that a node is selected (so that its neighboring edges will be corrupted)
%% p_edge_crpt: probability that an edge (among the neighboring edges of a fixed node) is corrupted
%% crpt_type: determines how we corrupt the relative permutation matrix given an edge to be corrupted
%  we have 4 options: 
% 'uniform': each corrupted relative permutations X_{ij} i.i.d uniformly drawn
% 'self-consistent': corrupted X_{ij} are rel. permutations of another set of ... 
%                    absolute permutations (different from ground truth).
%                    Namely X_{ij} = P_i^{adv} P_j^{adv}'
% 'local-biased': self-consistent corruption with additional sampling rejection procedure ...
%                 so that the expectation of rel. permutations deviates away from the ground truth
% 'local-adv': super malicious corruption that replaces the underlying absolute pemutation from 
%              P_i^{ground truth} to P_i^{adv} Namely X_{ij} = P_i^{adv} P_j^{ground truth}'. 
%              Addtional noise is added to the corruption, otherwise the recovery of ground truth is ill-posed.

%% Output:
%% AdjMat: n by n adjacency matrix of the graph G([n], E)
%% CrptMat: n by n adjacency matrix of the corrupted subgraph
%% Ind: the indices matrix of the graph G([n], E). Each row is a 2d vector (i,j), the index of an edge (i<j).
%%      the edge indices are sorted in row-major order, for example (1,2), (1,3),..., (2,4),(2,8),..
%% X: block matrix of relative permutations. Each [i,j]-th block is 
%%    the d by d relative permutation Xij between node i and j
%% X_orig: block matrix of ground truth relative permutations. Each [i,j]-th block is 
%%    the d by d true relative permutation Xij^* between node i and j

%% Reference
%% [1] Yunpeng Shi, Shaohan Li and Gilad Lerman. 
%% "Robust Multi-object Matching via Iterative Reweighting of the Graph Connection Laplacian" NeurIPS 2020.


function[AdjMat, CrptMat, Ind, X, X_orig] = Nonuniform_Topology(n,d,p, p_node_crpt,p_edge_crpt,crpt_type)

    G = rand(n,n) < p;
    G = tril(G,-1);
    % generate adjacency matrix
    AdjMat = G + G'; 
    [Ind_j, Ind_i] = find(G==1);
    Ind = [Ind_i, Ind_j];
    Ind_full = [Ind_j, Ind_i;Ind_i, Ind_j];
    m = length(Ind_i);

    P_orig = zeros(d,d,n);
    X_orig = zeros(n*d);
    for i = 1:n
        TempMat=eye(d);
        P_orig(:,:,i)=TempMat(randperm(d),:);
    end

    for k = 1:m
        i=Ind_i(k); j=Ind_j(k); 
        X_orig((d*i-(d-1)):(d*i), (d*j-(d-1)):(d*j))=P_orig(:,:,i)*(P_orig(:,:,j)');
    end
    X_orig = X_orig + X_orig';
    X = X_orig;

    node_crpt = randperm(n);
    n_node_crpt = floor(n*p_node_crpt);
    node_crpt = node_crpt(1:n_node_crpt);
    CrptMat = zeros(n,n);


    P_crpt = zeros(d,d,n);

    for i = 1:n
        TempMat=eye(d);
        P_crpt(:,:,i)=TempMat(randperm(d),:);
    end

    for i = node_crpt
        neighbor_cand = Ind_full(Ind_full(:,1)==i,2);
        neighbor_cand = reshape(neighbor_cand, 1, length(neighbor_cand));
        neighbor_crpt = randperm(length(neighbor_cand));
        n_neighbor = floor(p_edge_crpt * length(neighbor_cand));
        neighbor_crpt =    neighbor_crpt(1:n_neighbor);
        neighbor_crpt = neighbor_cand(neighbor_crpt);

        for j = neighbor_crpt 

            CrptMat(i,j) = 1;CrptMat(j,i) = 1;        
            TempMat=eye(d);
            P0=TempMat(randperm(d),:);

            if strcmp(crpt_type,'uniform')
                X((d*i-(d-1)):(d*i), (d*j-(d-1)):(d*j))= P0;    
                X((d*j-(d-1)):(d*j), (d*i-(d-1)):(d*i))= P0';
            end
            if strcmp(crpt_type,'self-consistent')
                X((d*i-(d-1)):(d*i), (d*j-(d-1)):(d*j))=P_crpt(:,:,i)*(P_crpt(:,:,j)');
                X((d*j-(d-1)):(d*j), (d*i-(d-1)):(d*i))=(P_crpt(:,:,i)*(P_crpt(:,:,j)'))';
            end
            if strcmp(crpt_type,'local-biased')
                if (sum((P_orig(:,:,i)'*P_crpt(:,:,i)).*((P_orig(:,:,j)'*P_crpt(:,:,j))),'all')>=2)
                    X((d*i-(d-1)):(d*i), (d*j-(d-1)):(d*j))= P0;  
                    X((d*j-(d-1)):(d*j), (d*i-(d-1)):(d*i))= P0';
                else
                    X((d*i-(d-1)):(d*i), (d*j-(d-1)):(d*j)) = P_crpt(:,:,i)*(P_crpt(:,:,j)');
                    X((d*j-(d-1)):(d*j), (d*i-(d-1)):(d*i)) = (P_crpt(:,:,i)*(P_crpt(:,:,j)'))';
                end
            end
            if strcmp(crpt_type,'local-adv')
                permind=randperm(d,3);
                P_temp=zeros(d,d);
                P_temp(setdiff(1:d,permind),setdiff(1:d,permind)) = eye(d-3);
                P_temp(permind,permind)=[0,1,0;0,0,1;1,0,0];
                X((d*i-(d-1)):(d*i), (d*j-(d-1)):(d*j))=P_temp*P_orig(:,:,j)';
                X((d*j-(d-1)):(d*j), (d*i-(d-1)):(d*i))=(P_temp*P_orig(:,:,j)')';
            end
        end     

    end


end
