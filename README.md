# Gromacs:
## Extract a frame:
1) gmx check -f traj.xtc
2) gmx trjconv -s topol.tpr -f traj.xtc -o frame223.gro -dump 5009

## Replace single ion with complex ions e.g. NO3:
1) using https://www.swissparam.ch/results.php?job=151171367970 and [O-][N+]([O-])=O, to generate forcefield
2) gmx editconf -f lig.pdb -o lig.gro
3) rearrange the force field
4) Remove the Cl- using parmed
5) gmx insert-molecules -f NoCl.gro -ci lig.gro -nmol 99 -try 200000

##Combing trajectory
gmx trjcat -f system_run001.xtc system_run002.xtc system_run003.xtc -o combined.xtc (-settime)
