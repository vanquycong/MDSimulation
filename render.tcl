update
update idle

set scaling 3

color Display Background white
display depthcue off
display shadows on
display ambientocclusion on
display projection Orthographic
axes location Off
color Name C silver
light 2 on
light 3 on

display resize [expr 850*$scaling] [expr 1000*$scaling]
display update ui

set mol [mol new barrel.psf]
mol addfile barrel.pdb $mol

#repman clear $mol
#repman add Membrane_$mol -molid $mol -sel "segname MEMB and noh and same residue as y>-10 and not residue 2071 2104 2074 2082 2075 2103" -coloring Name -style VDW -material AOChalky
#repman add Protein_$mol -molid $mol -sel "protein" -coloring "Structure" -style "NewCartoon resolution 100" -material "AOShiny"

molinfo $mol set center_matrix {{{1 0 0 17.283} {0 1 0 -10.917} {0 0 1 21.22} {0 0 0 1}}}
molinfo $mol set scale_matrix {{{0.0210307 0 0 0} {0 0.0210307 0 0} {0 0 0.0210307 0} {0 0 0 1}}}
molinfo $mol set rotate_matrix {{{-1 0 0 0} {0 0 -1 0} {0 -1 0 0} {0 0 0 1}}}
molinfo $mol set global_matrix {{{1 0 0 0.0308656} {0 1 0 -0.0857134} {0 0 1 0.14653} {0 0 0 1}}}

render TachyonLOptiXInternal barrel.png
