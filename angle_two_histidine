# ============================================================
# Calculate the angle between two histidine imidazole rings
# from a trajectory loaded in VMD.
#
# Output:
#   frame    theta(deg)
#
# theta = angle between the two imidazole ring-plane normals
# Range: 0 - 90 degrees
# ============================================================

# ---------- USER INPUT --------------------------------------

# Change these atom indices to your two imidazole rings.
# VMD atom indices start from 0.

set ring1 {1 2 3 277 278}
set ring2 {5 6 7 279 280}

# Output file
set outfile "ring_angle.dat"

# ------------------------------------------------------------

# Get the top molecule
set mol top

# Open output file
set fp [open $outfile w]

puts $fp "# frame    theta(deg)"

# Number of frames
set nframes [molinfo $mol get numframes]

# ------------------------------------------------------------
# Function: calculate plane normal from selected atoms
# ------------------------------------------------------------

proc plane_normal {sel} {

    # Center of geometry
    set center [measure center $sel weight none]

    # Move coordinates relative to center
    set coords {}

    foreach coord [$sel get {x y z}] {

        set x [expr {[lindex $coord 0] - [lindex $center 0]}]
        set y [expr {[lindex $coord 1] - [lindex $center 1]}]
        set z [expr {[lindex $coord 2] - [lindex $center 2]}]

        lappend coords [list $x $y $z]
    }

    # Use first three independent vectors to define plane
    set r1 [lindex $coords 0]
    set r2 [lindex $coords 1]
    set r3 [lindex $coords 2]

    # v1 = r2-r1
    set v1 [vecsub $r2 $r1]

    # v2 = r3-r1
    set v2 [vecsub $r3 $r1]

    # Normal = v1 x v2
    set normal [veccross $v1 $v2]

    # Normalize
    set normal [vecnorm $normal]

    return $normal
}

# ------------------------------------------------------------
# Loop over trajectory
# ------------------------------------------------------------

for {set frame 0} {$frame < $nframes} {incr frame} {

    # Go to frame
    animate goto $frame

    # Create atom selections
    set sel1 [atomselect $mol "index [join $ring1 { }]" frame $frame]
    set sel2 [atomselect $mol "index [join $ring2 { }]" frame $frame]

    # Calculate plane normals
    set n1 [plane_normal $sel1]
    set n2 [plane_normal $sel2]

    # Dot product
    set dot [vecdot $n1 $n2]

    # Because the direction of a plane normal is arbitrary,
    # use absolute value so theta is between 0 and 90 degrees.
    set dot [expr {abs($dot)}]

    # Protect against numerical roundoff
    if {$dot > 1.0} {
        set dot 1.0
    }

    # Angle in radians -> degrees
    set theta [expr {acos($dot) * 180.0 / acos(-1.0)}]

    # Write result
    puts $fp [format "%8d %12.5f" $frame $theta]

    # Delete selections
    $sel1 delete
    $sel2 delete
}

close $fp

puts "=============================================="
puts "Finished."
puts "Number of frames: $nframes"
puts "Output file: $outfile"
puts "=============================================="
