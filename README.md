# Robust Multi-object Matching via Iterative Reweighting of the Graph Connection Laplacian
This repo contains matlab files for implementing the method of the following paper
[Robust Multi-object Matching via Iterative Reweighting of the Graph Connection Laplacian](https://proceedings.neurips.cc/paper/2020/hash/ae06fbdc519bddaa88aa1b24bace4500-Abstract.html)

## Implementation
Put all the files in the same directory and run ``demo_IRGCL.m``. It calls functions ``Uniform_Toplogy.m'' or ''Nonuniform_Topology.m`` for generating the data that follow various corruption models. Then it calls function ``IRGCL.m`` to solve permutation synchronization. See comments in the begining of each functions for details.
