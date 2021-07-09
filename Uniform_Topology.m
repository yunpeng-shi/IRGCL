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
%  we have 2 options: 
% 'uniform': each corrupted relative permutations X_{ij} i.i.d uniformly drawn
% 'self-consistent': corrupted X_{ij} are rel. permutations of another set of ... 
%                    absolute permutations (different from ground truth).
%                    Namely X_{ij} = P_i^{adv} P_j^{adv}'


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

function[AdjMat, CrptMat, Ind, X, X_orig] = Uniform_Topology(n,d,p,q,crpt_type)

    G = rand(n,n) < p;
    G = tril(G,-1);
    % generate adjacency matrix
    AdjMat = G + G'; 
    [Ind_j, Ind_i] = find(G==1);
    Ind = [Ind_i, Ind_j];
    m = length(Ind_i);

    P_orig = zeros(d,d,n);
    X_orig = zeros(n*d);
    for i = 1:n
        TempMat=eye(d);
        P_orig(:,:,i)=TempMat(randperm(d),:);
    end

    for k = 1:m
        i=Ind_i(k); j=Ind_j(k); 
        X_orig((d*i-(d-1)):(d*i), (d*j-(d-1)):(d*j))=P_orig(:,:,i)*(P_orig(:,:,j))';
    end
    X = X_orig;
    X_orig = X_orig + X_orig';


    % indices of corrupted edges
    crptIndLog = (rand(1,m)<=q);
    crptInd=find(crptIndLog);
    CrptMat = zeros(n,n);
    if strcmp(crpt_type,'uniform')
        for k=crptInd
            i = Ind_i(k); j = Ind_j(k);
            CrptMat(i,j) = 1;CrptMat(j,i) = 1;    
            TempMat=eye(d);
            X((d*i-(d-1)):(d*i), (d*j-(d-1)):(d*j))=TempMat(randperm(d),:); 
        end
    end
    if strcmp(crpt_type,'self-consistent')
        P_crpt = zeros(d,d,n);
        for i = 1:n
            TempMat=eye(d);
            P_crpt(:,:,i)=TempMat(randperm(d),:);
        end
        for k=crptInd
        i = Ind_i(k); j = Ind_j(k);
        CrptMat(i,j) = 1;CrptMat(j,i) = 1;    
        X((d*i-(d-1)):(d*i), (d*j-(d-1)):(d*j))=P_crpt(:,:,i)*(P_crpt(:,:,j))';
        end
    end



    X = X+X';

end
