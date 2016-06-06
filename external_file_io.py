import glob
import utilities
import collections
import shutil
import numpy
import sys
from pdb_bfactor import *

import os
import subprocess
from subprocess import PIPE
import argparse

"""
FILE IO utilities
"""

def write_ptraj_align_file(ref, trjfile):
    ohandle=open('ptraj-align.in', 'w')
    ohandle.write('''
trajin {0}
reference {1}
rms reference out rmsd-ca.dat @CA,C,N,O
trajout aln-{0}'''.format(trjfile, ref))
    return


def write_rmsd_ptraj_file(ref, trjfile):
    ohandle=open('ptraj-rmsd.in', 'w')
    ohandle.write('''
trajin {0}
reference {1}
# must align by bb first
rmsd bb @CA reference out rmsd-aggregate-bb.dat perres perresout rmsd-bb-perresout.dat perresavg rmsd-bb-perresavg.dat
rmsd all !@H* reference out rmsd-aggregate-all.dat perres perresout rmsd-all-perresout.dat perresavg rmsd-all-perresavg.dat
atomicfluct out rmsf-all-perresavg.dat !@H* byres
atomicfluct out rmsf-bb-perresavg.dat @CA byres
'''.format(trjfile, ref))
    return


def write_strip_ptraj_file(ref, ambmask, outname):
    ohandle=open('ptraj-strip.in', 'w')
    ohandle.write('''
trajin  {0}
strip !{1}
trajout {2}-cluster-check.pdb'''.format(ref, ambmask, outname))
    return
    

def write_clustering_ptraj_file(ref, trjfile, ambmask, d, outname):
    ohandle=open('ptraj-cluster.in', 'w')
    ohandle.write('''
trajin  {1}
reference {0}
rms reference out rmsd-ca.dat @CA
cluster {2} mass epsilon {3} out {4}-cluster-d{3}_out averagelinkage gracecolor summary {4}-cluster-d{3}-summary_out info {4}-cluster-d{3}-Cluster_info repout {4}-cluster-d{3}-centroid repfmt pdb nofit
'''.format(ref, trjfile, ambmask, d, outname))
    return


def write_hbond_ptraj(trjfile, mask1, mask2, outname):
    ohandle=open('ptraj-hbonds.in', 'w')
    ohandle.write('''
trajin {0}
hbond H1 donormask {1}@N*,O* acceptormask {2}&@N*,O* nointramol out numhb-{3}-donor.dat  avgout avghb-{3}-donor.dat series 
hbond H2 acceptormask {1}&@N*,O* donormask {2}&@N*,O* nointramol out numhb-{3}-accept.dat avgout avghb-{3}-accept.dat series 
'''.format(trjfile, mask1, mask2, outname))
    return

def write_solvent_ptraj_file(ref, trjfile, mask, outname, occupancy):
#radial rad.dat 0.1 10.0 :T3P@O :346
#radial rad2.dat 0.1 10.0 :T3P@H1 :346 
#radial rad3.dat 0.1 10.0 :T3P@H2 :346 
    ohandle=open('ptraj-water.in', 'w')
    ohandle.write('''
trajin {0}
reference {1}
# Create average of solute to view with grid.
#center {2}
rms reference @C,CA,O,N
grid {3}_water.dx 50 0.5 50 0.5 50 0.5 :T3P@O pdb {3}_water{4}.pdb max {4}'''.format(trjfile, ref, mask, outname, occupancy))
    if mask!=None:
        ohandle.write('''
hbond H1 out {0}_hbond.dat {1}@N*,O* series avgout avgout.dat  printatomnum nointramol solventdonor :T3P@O  solventacceptor :T3P@O  solvout {0}_solvout.dat bridgeout {0}_bridgeout.dat 
run
runanalysis lifetime H1[solventhb] out {0}_lifetime.dat
'''.format(outname, mask))
    return


def write_load_bfactor_pml(datatype, ref_basename, maxval):
    maxval=int(numpy.ceil(maxval))
    ohandle=open('load_bfactor.pml', 'w')
    ohandle.write('''
load bfactor-%s-%s
spectrum b, blue_white_red, minimum=0, maximum=%s
show cartoon
hide lines
hide (h. and (e. c extend 1))
''' % (datatype, ref_basename, maxval))
    return

def write_all_hbonds_to_pml(files, outname, residue_mapper, ref):
    pymol_handle=open('%s_hbonds.pml' % outname, 'w')
    pymol_handle.write('load %s\n' % ref)
    pymol_handle.write('sel protein, polymer\n')
    pymol_handle.write('show cartoon, protein\n')
    pymol_handle.write('hide lines, protein\n')
    ohandle=open('%s_all_hbonds.dat' % outname, 'w')
    # need a damn reverse dict for this
    reverse_dict = {value: keypath for keypath, value in utilities.make_keypaths(residue_mapper)}
    for file in files:
        print "-----%s----" % file
        fhandle=open(file)
        for line in fhandle.readlines():
            if '#' in line:
                continue
            acceptor=line.split()[0] 
            donor=line.split()[1] 
            fraction=line.split()[4]
            percent=float(fraction)*100
            distance=line.split()[5]
            angle=line.split()[6]
            if fraction < 0.05:
                continue
            res1_chain, res1_name, orig_res1_num, atom1=utilities.parse_hbond(acceptor, reverse_dict)
            res2_chain, res2_name, orig_res2_num, atom2=utilities.parse_hbond(donor, reverse_dict)
            hbond='%s.%s%s@%s-%s.%s%s@%s' % (res1_chain, res1_name,orig_res1_num,atom1,res2_chain, res2_name, orig_res2_num,atom2)
            print hbond, percent
            ohandle.write('%s\t%s\n' % (hbond, percent))
            pymol_accept='%s/%s/%s' % (res1_chain, orig_res1_num, atom1)
            pymol_donor='%s/%s/%s' % (res2_chain, orig_res2_num, atom2)
            pymol_handle.write('show sticks, chain %s and resi %s\n' % (res1_chain, orig_res1_num))
            pymol_handle.write('show sticks, chain %s and resi %s\n' % (res2_chain, orig_res2_num))
            hcolor=utilities.percent_score(percent)
            if hcolor=='red':
                pymol_handle.write('dist %s_strong_hbonds, %s, %s\n' % (outname, pymol_donor, pymol_accept))
                pymol_handle.write('color %s, %s_strong_hbonds\n' % (hcolor, outname))
            if hcolor=='pink':
                pymol_handle.write('dist %s_medium_hbonds, %s, %s\n' % (outname, pymol_donor, pymol_accept))
                pymol_handle.write('color %s, %s_medium_hbonds\n' % (hcolor, outname))
            if hcolor=='yellow':
                pymol_handle.write('dist %s_weak_hbonds, %s, %s\n' % (outname, pymol_donor, pymol_accept))
                pymol_handle.write('color %s, %s_weak_hbonds\n' % (hcolor, outname))
    pymol_handle.write('hide labels, %s_weak_hbonds\n' % outname)
    pymol_handle.write('hide labels, %s_medium_hbonds\n' % outname)
    pymol_handle.write('hide labels, %s_strong_hbonds\n' % outname)
    pymol_handle.close()
    ohandle.close()
    return 

def write_bridged_residues_to_pml(file, amberdata, total_frames, outname, reverse_dict, ref):
    pymol_handle=open('%s_bridged.pml' % outname, 'w')
    pymol_handle.write('load %s\n' % ref)
    pymol_handle.write('sel protein, polymer\n')
    pymol_handle.write('show cartoon, protein\n')
    pymol_handle.write('hide lines, protein\n')
    fhandle=open(file)
    for line in fhandle.readlines():
        if '#' in line:
            pass
        else:
            frames=int(line.split()[5])
            fraction=float(frames)/total_frames
            percent=fraction*100
            if percent < 5:
                break
            res1_num=line.split()[2].split(':')[0]
            res1_name=line.split()[2].split(':')[1]
            res2_num=line.split()[3].split(':')[0]
            res2_name=line.split()[3].split(':')[1]
            res1_chain=reverse_dict[res1_num][0]
            orig_res1_num=reverse_dict[res1_num][1]
            res2_chain=reverse_dict[res2_num][0]
            orig_res2_num=reverse_dict[res2_num][1]
            # determine the atoms from other hbond data, le sigh
            location=numpy.abs(numpy.array(amberdata[res1_num]['percent'])-percent).argmin()
            atom1=amberdata[res1_num]['atom'][location]
            atom1=atom1.split('_')[0]
            location=numpy.abs(numpy.array(amberdata[res2_num]['percent'])-percent).argmin()
            atom2=amberdata[res2_num]['atom'][location]
            atom2=atom2.split('_')[0]
            #print 'bridge persists for %0.2f %% of simulation' % percent
            pymol_res1='%s/%s/%s' % (res1_chain, orig_res1_num, atom1)
            pymol_res2='%s/%s/%s' % (res2_chain, orig_res2_num, atom2)
            pymol_handle.write('show sticks, chain %s and resi %s\n' % (res1_chain, orig_res1_num))
            pymol_handle.write('show sticks, chain %s and resi %s\n' % (res2_chain, orig_res2_num))
            hcolor=utilities.percent_score(percent)
            if hcolor=='red':
                pymol_handle.write('dist %s_strong_bridge, %s, %s\n' % (outname, pymol_res1, pymol_res2))
                pymol_handle.write('color %s, %s_strong_bridge\n' % (hcolor, outname))
            if hcolor=='pink':
                pymol_handle.write('dist %s_medium_bridge, %s, %s\n' % (outname, pymol_res1, pymol_res2))
                pymol_handle.write('color %s, %s_medium_bridge\n' % (hcolor, outname))
            if hcolor=='yellow':
                pymol_handle.write('dist %s_weak_bridge, %s, %s\n' % (outname, pymol_res1, pymol_res2))
                pymol_handle.write('color %s, %s_weak_bridge\n' % (hcolor, outname))
    pymol_handle.write('hide labels, %s_weak_bridge\n' % outname)
    pymol_handle.write('hide labels, %s_medium_bridge\n' % outname)
    pymol_handle.write('hide labels, %s_strong_bridge\n' % outname)
    fhandle.close()
    pymol_handle.close()
    return 

