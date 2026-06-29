# Code
## Extract a frame:
1) gmx check -f traj.xtc
2) gmx trjconv -s topol.tpr -f traj.xtc -o frame223.gro -dump 5009
