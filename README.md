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

## Example execution:
`./protein_interface.py -i 4hhb.pdb -c1 A -c2 B -t 7

The PDB input file given is: 4hhb.pdb
First chain (Chain1) specified is: A
Second chain (Chain2) specified is: B
The threshold distance specifies is: 7


A:GLU(30) interacts with B:PRO(124)
A:ARG(31) interacts with B:PRO(124)
A:LEU(34) interacts with B:PRO(124)
A:LEU(34) interacts with B:PRO(125)
A:LEU(34) interacts with B:ALA(128)
A:SER(35) interacts with B:GLN(127)
A:SER(35) interacts with B:ALA(128)
A:SER(35) interacts with B:GLN(131)
A:LEU(106) interacts with B:CYS(112)
A:VAL(107) interacts with B:CYS(112)
A:VAL(107) interacts with B:ALA(115)
A:ALA(110) interacts with B:CYS(112)
A:ALA(110) interacts with B:ALA(115)
A:ALA(110) interacts with B:HIS(116)
A:ALA(111) interacts with B:ALA(115)
A:ALA(111) interacts with B:HIS(116)
A:ALA(111) interacts with B:GLY(119)
A:ALA(111) interacts with B:LYS(120)
A:HIS(112) interacts with B:GLY(119)
A:PRO(114) interacts with B:HIS(116)
A:PRO(119) interacts with B:ARG(30)
A:PRO(119) interacts with B:VAL(33)
A:PRO(119) interacts with B:PRO(51)
A:ALA(120) interacts with B:PRO(51)
A:HIS(122) interacts with B:VAL(34)
A:ALA(123) interacts with B:VAL(33)
A:ALA(123) interacts with B:VAL(34)


A chain
13/14 of the interface amino acids lying on alpha helices.
0/14 of the interface amino acids lying on beta sheets.

B chain
12/14 of the interface amino acids lying on alpha helices.
0/14 of the interface amino acids lying on beta sheets.

A chain
GLU: closest ARG at distance 1
ARG: closest LEU at distance 3
LEU: closest SER at distance 1
SER: closest LEU at distance 71
LEU: closest VAL at distance 1
VAL: closest ALA at distance 3
ALA: closest ALA at distance 1
ALA: closest HIS at distance 1
HIS: closest PRO at distance 2
PRO: closest PRO at distance 5
PRO: closest ALA at distance 1
ALA: closest HIS at distance 2
HIS: closest ALA at distance 1

B chain
ARG: closest VAL at distance 3
VAL: closest VAL at distance 1
VAL: closest PRO at distance 17
PRO: closest CYS at distance 61
CYS: closest ALA at distance 3
ALA: closest HIS at distance 1
HIS: closest GLY at distance 3
GLY: closest LYS at distance 1
LYS: closest PRO at distance 4
PRO: closest PRO at distance 1
PRO: closest GLN at distance 2
GLN: closest ALA at distance 1
ALA: closest GLN at distance 3`
