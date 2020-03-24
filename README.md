# Protein_Interface
Python script for determining interface atoms between two chains of a protein using Euclidean distance.<br/>
Many biological functions involve the formation of protein-protein complexes. A complex is a structure consisting of two or more interacting proteins. The interface of interacting proteins is defined as the set of amino acids of the two proteins within a certain threshold distance. The characterization of protein-protein  interfaces is important because it will enable the prediction of protein interactions providing insight into their function.

## Input arguments:
- -i: Input PDB file
- -c1: Chain1 of PDB file
- -c2: Chain2 of PDB file
- t: Distance threshold

## Execution of script:
`./protein_interface.py -i <input PDB file> -c1 <Chain 1> -c2 <Chain 2> -t <distance threshold>`
