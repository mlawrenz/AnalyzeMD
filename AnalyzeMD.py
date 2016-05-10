from schrodinger.trajectory.desmondsimulation import create_simulation
from pdb_bfactor import *
import sys
import schrodinger.trajectory.analysis as analysis
from schrodinger.structure import write_ct, write_ct_pdb
from schrodinger.structutils.analyze import evaluate_asl
from schrodinger.application.desmond.util import get_indices, parse_slice
import os
import subprocess
from subprocess import PIPE
import argparse

"""
Analyze MD
===========

** To be used on Schrodinger Desmond results  **

Program requires as input (see args with -h flag):
*
*

Program requires AMBERTools and will produce DCD files
This file can be loaded into VMD:
* you can run the command: vmd -pdb reference.pdb -dcd file.dcd

"""

def run_linux_process(command):
    p=subprocess.Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    p.wait()
    output, err=p.communicate()
    return output, err

def write_rmsd_ptraj_file(ref, trjfile):
    ohandle=open('ptraj-rmsd.in', 'w')
    ohandle.write('''
trajin {0}
reference {1}
# must align by bb first
rmsd bb @CA reference out rmsd-bb.dat perres perresout rmsd-bb-perresout.dat perresavg rmsd-bb-perresavg.dat
rmsd all !@H* reference out rmsd-all.dat perres perresout rmsd-all-perresout.dat perresavg rmsd-all-perresavg.dat
atomicfluct out rmsf-all.dat !@H* byres
atomicfluct out rmsf-bb.dat @CA byres
'''.format(trjfile, ref))
    return


def write_clustering_ptraj_file(ref, trjfile, ambmask, d):
    ohandle=open('ptraj-cluster.in', 'w')
    ohandle.write('''
trajin  {0}
reference {1}
rms reference out rmsd-ca.dat @CA
cluster {2} mass epsilon {3} out cluster-d{3}_out averagelinkage gracecolor summary cluster-d{3}-summary_out info cluster-d{3}-Cluster_info repout cluster-d{3}-centroid repfmt pdb
'''.format(ref, trjfile, ambmask, d))
    return


def write_load_bfactor_pml(ref_basename):
    ohandle=open('load_bfactor.pml', 'w')
    ohandle.write('''
load bfactor-%s
spectrum b, blue_white_red, minimum=0, maximum=3
show cartoon
hide lines
hide (h. and (e. c extend 1))
''' % ref_basename)
    return

def make_analysis_folder(cwd, name):
    if not os.path.exists('%s/%s-analysis' % (cwd, name)):
        os.mkdir('%s/%s-analysis' % (cwd, name))
    os.chdir('%s/%s-analysis' % (cwd, name))
    return

def rmsd_and_rmsf(cwd, ref, trjfile, datatype):
    ref=os.path.abspath(ref)
    ref_basename=os.path.basename(ref)
    trjfile=os.path.abspath(trjfile)
    #make sure have absolute path since working in analysis folder now
    make_analysis_folder(cwd, datatype)
    write_rmsd_ptraj_file(ref, trjfile)
    command='%s/bin/cpptraj %s ptraj-rmsd.in' %  (os.environ['AMBERHOME'], ref)
    output, err=run_linux_process(command)
    if 'rror' in err:
        print "ERROR IN CPPTRAJ RMSD CALCULATION"
        print err
        sys.exit()
    f = open(ref,'r')
    pdb_data = f.readlines()
    f.close()
    if datatype=='rmsd-all':
        datafile='rmsd-all-perresavg.dat'
    if datatype=='rmsd-bb':
        datafile='rmsd-bb-perresavg.dat'
    if datatype=='rmsf-all':
        datafile='rmsf-all.dat'
    if datatype=='rmsf-bb':
        datafile='rmsf-bb.dat'
    out=map_file_by_index(datafile, pdb_data)
    print "Warning: Mapping data by residue index from datafile to ref pdbfile"
    print "Warning: Make sure atoms correpond!"
    outfilename='bfactor-%s' % ref_basename
    outfile=open(outfilename, 'w')
    for line in out:
        outfile.write(line)
    outfile.close()
    write_load_bfactor_pml(ref)

def hbond_analysis():
    return

def clustering(ref, trjfile, cluster, d=1.0):
    ambmask=':%s' % cluster
    write_clustering_ptraj_file(ref, trjfile, ambmask, d)
    command='%s/bin/cpptraj %s ptraj-cluster.in' %  (os.environ['AMBERHOME'], ref)
    output, err=run_linux_process(command)


def main(args):
    cwd=os.getcwd()
    if not os.environ['AMBERHOME']:
        print "AMBERHOME IS NOT SET"
        sys.exit()
    if args.analysis=='rmsd_calc':
        print "RUNNING RMSD/RMSF ANALYSIS"
        rmsd_and_rmsf(cwd, args.reffile, args.trjfile, args.rmsdtype)
    if args.analysis=='clustering':
        print "RUNNING CLUSTERING ON RESIDUES %s" % args.clustermask
        clustering(cwd, args.reffile, args.trjfile, args.clustermask, args.distance)
    print "done"
    


def parse_commandline():
    parser = argparse.ArgumentParser(description='''
Run chosen analysis on MD DCD
trajectory file. Choose rmsd or clustering. See additional options for each
analysis choice by running: 
$SCHRODINGER/run /home/mlawrenz/AnalyzeMD/AnalyzeMD.py clustering -h''')


    #parser.add_argument("analysis", choices=["rmsd", "clustering", "hbond"])
    parser.add_argument('-r', '--reffile', dest='reffile',
                      help='reference structure PDB file, for RMSD and visualization. MUST MATCH TRAJECTORY')
    parser.add_argument('-t', '--trjfile', dest='trjfile',
                      help='DCD trj for trajectory analysis')
    subparsers= parser.add_subparsers(help='analysis suboption choices', dest='analysis')

    r_parser=subparsers.add_parser("rmsd_calc")
    r_parser.add_argument('rmsdtype',choices=['rmsd-all', 'rmsd-bb', 'rmsf-all', 'rmsf-bb'], help='type of rmsd or rmsf calculation, using all heavy atoms or just backbone CA atoms')
    c_parser=subparsers.add_parser("clustering")
    c_parser.add_argument('-c', '--cluster', dest='clustermask',
                      help='residues to use for clustering, separated by commas and dashes i.e. 45-55,68,128,155')
    c_parser.add_argument('-d', '--distance', dest='distance', default=2.0,
                      help='clustering RMSD cutoff to define clusters. recommend 2 for most small changes.')
    args = parser.parse_args()
    return args

#run the function main if namespace is main
if __name__ == "__main__":
    args = parse_commandline()
    main(args)



