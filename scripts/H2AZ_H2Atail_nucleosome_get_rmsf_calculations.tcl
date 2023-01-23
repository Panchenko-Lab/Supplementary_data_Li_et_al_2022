
mol addfile H2AZ_H2Atail_nucleosome_initial.gro 
mol addfile H2AZ_H2Atail_nucleosome_overal_2us_long_final.pdb.trr step 1 
mol addfile H2AZ_H2Atail_nucleosome_average.gro 

set nframes [expr  [molinfo top get numframes] - 1 ]
set all [atomselect top "all"]
set rmsf [measure rmsf $all first 2001 last 20000]

$all frame $nframes
$all set beta $rmsf

set chains [list A E B F C G D H I J]
set chains_names [list H3_1 H3_2 H4_1 H4_2 H2A_1 H2A_2 H2B_1 H2B_2 DNA_I DNA_J]

set outfile [open "H2AZ_H2Atail_nucleosome_get_rmsf_calculations.dat"  w]

set ts "Resid"
foreach chain  $chains_names {
append ts "\t$chain"
}
puts $outfile $ts

for  { set i -73 } { $i<=150 } { incr i } {
set ds "$i"
foreach seg  $chains {

set sel [atomselect top "chain $seg and name CA P and resid '$i' "]
if {[$sel num] == 0} { set beta 0
} else {
set  beta [$sel get beta]
}

append ds "\t$beta"

}
puts $outfile $ds
}

close $outfile 


set ts "index"
foreach chain  $chains_names {
append ts "\t$chain"
}
puts $outfile $ts

for  { set i -73 } { $i<=150 } { incr i } {
set ds "$i"
foreach seg  $chains {
set sel [atomselect top "chain $seg and (not backbone) and resid '$i' and noh"]
set sel [atomselect top "chain $seg and (not backbone) and (not name P O1P O2P O3' O5' C5' C1' C2' C3' C4' O4') and resid '$i' and noh"]
if {[$sel num] == 0} { set mean 0
} else {
set  beta [$sel get beta]
set bl [llength $beta]
set mean 0
for {set j 0} { $j<$bl } { incr j } {
set mean [expr $mean + [lindex $beta $j]]
}
set mean [expr $mean / $bl]
}


append ds "\t$mean"

}
puts $outfile $ds
}

close $outfile 


set nframes [expr  [molinfo top get numframes] - 1 ]

puts -nonewline $outfile "\n"
for { set i 1 } { $i<=$nframes } { incr i } {
set time [expr 0.1 * $i]
puts -nonewline $outfile [format "%.3f" "$time"]
puts  [format "Time %.3f" "$time"]
for { set r1 -73 } { $r1<=73 } { incr r1 } {
set sel1 [atomselect top "protein and noh" frame $i]
set r2 [expr $r1 * (-1)]
set sel2 [atomselect top "((chainI and resid '$r1') or (chainJ and resid '$r2')) and noh" frame $i]
set contacts [lindex [measure contacts  3.0 $sel1 $sel2] 1]
}

foreach ind [$h3 get index] {

set sel0 [atomselect top "index $ind" frame 0]
set seli [atomselect top "index $ind" frame $i]

set c0 [lindex [$sel0 get {x y z}] 0]
set c1 [lindex [$seli get {x y z}] 0]

set d [veclength [vecsub $c0 $c1]]
puts $outfile "$d"

$sel0 delete
$seli delete
}
exit


