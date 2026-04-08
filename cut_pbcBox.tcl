 ================================
# USER PARAMETERS
# ================================

set resid_selection "resname HSE HSP"   ;# change if needed (e.g. "resid 128 199")
set box_half_length 8.0             ;# in nm (7.5 nm = 15 Å box)

set output_file "qm_cluster.gro"

# ================================
# LOAD pbctools
# ================================
package require pbctools

# ================================
# FIX PERIODICITY
# ================================
pbc unwrap -all
pbc wrap -all -compound res -center com

# ================================
# SELECT HISTIDINES AND GET CENTER
# ================================
set sel [atomselect top $resid_selection]
set center [measure center $sel weight mass]

puts "Center of selection: $center"

# ================================
# DEFINE CUBE AROUND CENTER
# ================================
set xmin [expr [lindex $center 0] - $box_half_length]
set xmax [expr [lindex $center 0] + $box_half_length]

set ymin [expr [lindex $center 1] - $box_half_length]
set ymax [expr [lindex $center 1] + $box_half_length]

set zmin [expr [lindex $center 2] - $box_half_length]
set zmax [expr [lindex $center 2] + $box_half_length]

# ================================
# SELECT ATOMS INSIDE CUBE
# KEEP WHOLE RESIDUES (important for water!)
# ================================
set cube_sel [atomselect top \
    "same residue as (x > $xmin and x < $xmax and y > $ymin and y < $ymax and z > $zmin and z < $zmax)"]

puts "Atoms selected: [$cube_sel num]"

# ================================
# WRITE OUTPUT
# ================================
$cube_sel writegro $output_file

puts "Done! Output written to $output_file"

puts "Done! Output written to $output_file"

# ================================
# CLEAN UP
# ================================
$sel delete
$cube_sel delete
