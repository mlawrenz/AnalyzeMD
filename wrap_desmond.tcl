mol new XXX-out.cms type mae first 0 last -1 step 1 waitfor all
mol addfile XXX_trj/clickme.dtr type dtr first 0 last -1 step 5 waitfor all
package require pbctools 
set solute [ atomselect top "protein or resname UNK" ]
set cell [ pbc get -now ]
pbc set $cell
pbc wrap -centersel "resname YYY" -center com -compound residue -all 
#pbc wrap -centersel "protein or resname UNK" -center com -compound residue -all 
#pbc wrap -centersel "protein" -center com -compound residue -all 
#pbc wrap -compound fragment -center com -all
set n [molinfo top get numframes]
for {set t 0} {$t < $n} {incr t} {
set all [atomselect top "all" frame $t]
set ref [atomselect top "backbone" frame 0]
set bb [atomselect top "backbone" frame $t]
set trans_mat [measure fit $bb $ref]
$all move $trans_mat
}
animate write dcd aln-wrap-XXX.dcd sel $all waitfor all 
animate goto start
animate write pdb aln-wrap-XXX.pdb sel $all beg 0 end 1
exit

