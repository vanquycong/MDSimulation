proc enabletrace {} {
    global vmd_frame
    trace variable vmd_frame([molinfo top]) w drawcounter
}

proc disabletrace {} {
    global vmd_frame
    trace vdelete vmd_frame([molinfo top]) w drawcounter
}

proc drawcounter {name element op} {
    global vmd_frame

    set mol [molinfo top]

    # Delete only the previous time label and not all the drawn objects
    if {[info exists ::time_label_id]} {
        draw delete $::time_label_id
    }

    draw color black

    set psperframe 24
    set psoffset 0

    set time [format "%d fs" \
        [expr ($vmd_frame($mol) * $psperframe) + $psoffset]]

    # Store the ID of the newly drawn text
    set ::time_label_id [draw text {46 44 46} "$time" \
        thickness 3 size 1]
}

enabletrace
