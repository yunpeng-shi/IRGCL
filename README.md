# Robust Multi-object Matching via Iterative Reweighting of the Graph Connection Laplacian
This repo contains matlab files for implementing the permutation synchronization solver of the following paper
[Robust Multi-object Matching via Iterative Reweighting of the Graph Connection Laplacian](https://proceedings.neurips.cc/paper/2020/hash/ae06fbdc519bddaa88aa1b24bace4500-Abstract.html). Yunpeng Shi, Shaohan Li, Gilad Lerman. NeurIPS 2020.

## Implementation
Put all the files in the same directory and run ``demo_IRGCL.m``. It calls functions ``Uniform_Topology.m`` or ``Nonuniform_Topology.m`` for generating the data that follow various corruption models. Then it calls function ``IRGCL.m`` to solve permutation synchronization. See comments in the begining of each function for details.

## A Variety of Corruption Models
We provide 6 different corruption models. 4 for nonuniform topology and 2 for uniform toplogy. Uniform/Nonuniform toplogy refers to whether the corrupted subgraph is Erdos Renyi or not. In other words, the choice of Uniform/Nonuniform toplogy decides how to select edges for corruption. In ``Uniform_Topology.m``, two nodes are connected with probability ``p``. Then edges are independently drawn with probability ``q`` for corruption. In ``Nonuniform_Topology.m``, two nodes are connected with probability ``p``. Then with probability ``p_node_crpt`` a node is selected so that its neighboring edges will be corrupted. Next, for each selected node, with probability ``p_edge_crpt`` an edge (among the neighboring edges of the selected node) is corrupted. 

The argument ``crpt_type`` in the two functions determines how the corrupted relative permutations are generated for those selected edges. In ``Uniform_Topology.m``, there are 2 options of ``crpt_type``: ``uniform`` and ``self-consistent``.
In ``Nonuniform_Topology.m``, there are the following 4 options of ``crpt_type``.

``uniform``: The corrupted relative permutations <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{X_{ij}}"> are i.i.d follows uniform distribution over the space of permutation group.

``self-consistent``: The corrupted <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{X_{ij}}"> are relative permutations of another set of absolute permutations. Namely <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{X_{ij} = P_i^{crpt} P_j^{crpt}'}"> where those absolute permutations are different from the ground truth and are i.i.d drawn from the uniform distribution in the space of permutation group.

``local-biased``: Self-consistent corruption with additional sampling rejection procedure so that the expectation of relative permutations deviates away from the ground truth.

``local-adv``: Extremely malicious corruption that replaces the underlying absolute pemutations from ground truth <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{P_i^*}"> to <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{P_i^{crpt}}">. Namely <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{X_{ij} = P_i^{crpt} P_j^{* }'}">. Additional high noise is added to the corruption, otherwise the recovery of the ground truth can be ill-posed.

