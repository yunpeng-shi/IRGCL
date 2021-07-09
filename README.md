# Robust Multi-object Matching via Iterative Reweighting of the Graph Connection Laplacian
This repo contains matlab files for implementing the method of the following paper
[Robust Multi-object Matching via Iterative Reweighting of the Graph Connection Laplacian](https://proceedings.neurips.cc/paper/2020/hash/ae06fbdc519bddaa88aa1b24bace4500-Abstract.html)

## Implementation
Put all the files in the same directory and run ``demo_IRGCL.m``. It calls functions ``Uniform_Toplogy.m`` or ``Nonuniform_Topology.m`` for generating the data that follow various corruption models. Then it calls function ``IRGCL.m`` to solve permutation synchronization. See comments in the begining of each function for details.

## Various Corruption Models
We provide 6 different corruption models. 4 for nonuniform topology and 2 for uniform toplogy. Uniform/Nonuniform toplogy refers to whether the corrupted subgraph is Erdos Renyi or not. In other words, the choice of Uniform/Nonuniform toplogy decides how to select edges for corruption. The argument ``crpt_type`` in the two functions decide how the corrupted relative permutations are generated for those selected edges.

There are 4 options of ``crpt_type`` argument for ``Nonuniform_Topology.m``.
``uniform``: each corrupted relative permutations $X_{ij}$ i.i.d uniformly drawn
``self-consistent``: corrupted $X_{ij}$ are rel. permutations of another set of absolute permutations (different from ground truth). Namely $X_{ij} = P_i^{crpt} P_j^{crpt}'$
``local-biased``: self-consistent corruption with additional sampling rejection procedure so that the expectation of rel. permutations deviates away from the ground truth.
``local-adv``: super malicious corruption that replaces the underlying absolute pemutation from ground truth $P_i^* $ to $P_i^{crpt}$. Namely $X_{ij} = P_i^{crpt} P_j^{* }'$. Addtional noise is added to the corruption, otherwise the recovery of ground truth is ill-posed.

