# Gromacs:
## Commands
### Extract a frame:
1) gmx check -f traj.xtc
2) gmx trjconv -s topol.tpr -f traj.xtc -o frame223.gro -dump 5009

### Combing trajectory
gmx trjcat -f system_run001.xtc system_run002.xtc system_run003.xtc -o combined.xtc (-settime)

## Neutralization
### Way1
1) Manually Creating "ions.mdp"
``
integrator = steep
emtol = 1000.0
emstep = 0.01
nsteps = 1 ``

3) gmx_mpi_d grompp -f ions.mdp -c system.gro -p system.top -o ions.tpr
4) Manually Changing "WAT" to "SOL" in [ molecules] of "system.top"
5) gmx_mpi_d genion -s ions.tpr -o system.gro -p system.top -pname SOD -nname CLA -neutral (-conc 0.15) 
6) Manually Changing "SOL" to "WAT" in [ molecules] of "system.top", because "SOL" is not in .gro file

### Way2
1) Manually Creating "ions.mdp"
2) gmx_mpi_d grompp -f ions.mdp -c system.gro -p system.top -o ions.tpr
3) CREATING Index Files by coping the SOL group and change it to WAT.
   gmx_mpi_d make_ndx -f ions.tpr
   19
   name 19 WAT
4) gmx_mpi_d genion -s ions.tpr -o system.gro -p system.top -pname SOD -nname CLA -neutral -conc 0.15 -n index.ndp
5) Choose: "19"

## Replace single ion with complex ions e.g. NO3:
1) using https://www.swissparam.ch/results.php?job=151171367970 and [O-][N+]([O-])=O, to generate forcefield
2) gmx editconf -f lig.pdb -o lig.gro
3) rearrange the force field
4) Remove the Cl- using parmed
5) gmx insert-molecules -f NoCl.gro -ci lig.gro -nmol 99 -try 200000
