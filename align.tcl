
mol load pdb aln-wrap-Cdc34A_acid_SS_50ns.pdb
mol addfile equil-combine-solvent-Cdc34A_acid_SS.dcd type dcd waitfor all
set n [molinfo top get numframes]
for {set t 0} {$t < $n} {incr t} {
set all [atomselect top "all" frame $t]
set ref [atomselect top "backbone" frame 0]
set bb [atomselect top "backbone" frame $t]
set trans_mat [measure fit $bb $ref]
$all move $trans_mat
}
set sel [atomselect top all]
set origin [measure center $sel]
set fd [open "origin.dat" w]
puts $fd $origin
close $fd

animate write dcd morganaligned-equil-combine-solvent-Cdc34A_acid_SS.dcd  sel $all waitfor all 
animate goto start
animate write pdb morganaligned-aln-wrap-Cdc34A_acid_SS_50ns.pdb sel $all beg 0 end 1
exit