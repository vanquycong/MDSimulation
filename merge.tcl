package require topotools

##################################### Input residue number##########
while {1} {
    puts "Enter the residue number for merge.tlc (e.g., 87):"
    flush stdout
    set resid [gets stdin]
    
    # Check if the input is a valid integer
    if {[string is integer -strict $resid]} {
        break ;# Exit the loop if the input is an integer
    } else {
        puts "Invalid input. Please enter a valid residue number."
    }
}




##################################################
set prot [mol new premerge.psf]
mol addfile premerge.pdb
set sel_prot [atomselect $prot "resid $resid"]

set rest [mol new premut_rest.psf]
mol addfile premut_rest.pdb

set final [::TopoTools::mergemols "$prot $rest"]
set sel [atomselect $final "all"]
set sel_system [atomselect $final "protein and resid $resid"]

$sel_system set beta [$sel_prot get beta] 

$sel writepsf fep.psf
$sel writepdb fep.pdb

quit

