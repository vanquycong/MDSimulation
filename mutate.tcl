package require mutator
# Prompt the user to input the amino acid
# Prompt the user to input the amino acid
while {1} {
    puts "Enter the amino acid (e.g. ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HSD, HSE, HSP, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL):"
    flush stdout
    set amino [gets stdin]
    
    # Check if the input is a valid string
    if {[string is string $amino]} {
        break ;# Exit the loop if the input is a string
    } else {
        puts "Invalid input. Please enter a valid amino acid."
    }
}


# Prompt the user to input the residue number
while {1} {
    puts "Enter the residue number (e.g., 87):"
    flush stdout
    set resid [gets stdin]
    
    # Check if the input is a valid integer
    if {[string is integer -strict $resid]} {
        break ;# Exit the loop if the input is an integer
    } else {
        puts "Invalid input. Please enter a valid residue number."
    }
}

# Load PSF and coordinate files
mol new system.psf
mol addfile system_run100.coor

# Create selections for protein and not protein
set sel_protein [atomselect top protein]
$sel_protein writepsf premut.psf
$sel_protein writepdb premut.pdb

set sel_not_protein [atomselect top "not protein"]
$sel_not_protein writepsf premut_rest.psf
$sel_not_protein writepdb premut_rest.pdb

# Perform mutation using Mutator plugin
mutator -psf premut.psf -pdb premut.pdb -o mut -resid $resid -mut $amino -FEP protein

# Check if the destination file exists
if {[file exists "premerge.pdb"]} {
    # Delete the existing file
    file delete -force "premerge.pdb"
}

if {[file exists "premerge.psf"]} {
    # Delete the existing file
    file delete -force "premerge.psf"
}

file rename "protein.fep" "premerge.pdb"
file rename "protein.fep.psf" "premerge.psf"

source "merge.tcl"


