import glob
import external_file_io
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
    make_analysis_folder(cwd, 'rmsd')
    external_file_io.write_rmsd_ptraj_file(ref, trjfile)
    command='%s/bin/cpptraj %s ptraj-rmsd.in' %  (os.environ['AMBERHOME'], ref)
    output, err=run_linux_process(command)
    numpy.savetxt('ptraj-rmsd.out', output.split('\n'), fmt='%s')
    if 'rror' in err:
        numpy.savetxt('ptraj-rmsd.err', err.split('\n'), fmt='%s')
        print "ERROR IN CPPTRAJ RMSD CALCULATION"
        print "CHECK ptraj-rmsd.err"
        print err
        sys.exit()
    f = open(ref,'r')
    pdb_data = f.readlines()
    f.close()
    datafile='%s-perresavg.dat' % datatype
    if datatype=='rmsf-bb':
        bb=check_pdb_file(pdb_data)
        maxval, out=map_file_by_index(datafile, pdb_data, bb)
    else:
        maxval, out=map_file_by_index(datafile, pdb_data)
    outfilename='bfactor-%s-%s' % (datatype, ref_basename)
    outfile=open(outfilename, 'w')
    for line in out:
        outfile.write(line)
    outfile.close()
    external_file_io.write_load_bfactor_pml(datatype, ref_basename, maxval)
    read_handle=open('README', 'w')
    read_handle.write('''
rmsd-aggregate-all.dat  # RMSD of all atoms in selection for each time point
rmsd-all-perresavg.dat  # RMSD average for each residue over the simulation
rmsd-all-perresout.dat  # RMSD for each residue for each time point
rmsf-all-perresavg.dat  # RMSF average for each residue over the simulation
bfactor-%s-%s           # PDB file with B factors filled with specified datatype
''' % (datatype, ref_basename))
    read_handle.close()
    return

def sovlent_calc(cwd, outname, ref, trjfile, selection=None, occupancy=0.8):
    ref=os.path.abspath(ref)
    ref_basename=os.path.basename(ref)
    trjfile=os.path.abspath(trjfile)
    residue_mapper=utilities.map_residues(ref)
    reverse_dict = {value: keypath for keypath, value in utilities.make_keypaths(residue_mapper)}
    #make sure have absolute path since working in analysis folder now
    make_analysis_folder(cwd, 'solvent')
    # test that the residues you cluster on are the right ones
    if selection!=None:
        new_ambmask=utilities.check_cluster_pdbfile(ref, selection, outname)
        external_file_io.write_solvent_ptraj_file(ref, trjfile, new_ambmask, outname, occupancy)
    else:
        external_file_io.write_solvent_ptraj_file(ref, trjfile, selection, outname, occupancy)
    command='%s/bin/cpptraj %s ptraj-water.in' %  (os.environ['AMBERHOME'], ref)
    output, err=run_linux_process(command)
    if 'rror' in err:
        numpy.savetxt('ptraj-water.err', err.split('\n'), fmt='%s')
        print "ERROR IN CPPTRAJ RMSD CALCULATION"
        print "CHECK ptraj-water.err"
        print err
        sys.exit()
    print "PDB file with {1} occupancy are in {0}_water{1}.pdb".format(outname, occupancy)

    if selection!=None:
        file='%s_solvout.dat' % outname
        fhandle=open(file)
        for line in fhandle.readlines():
            if '#' in line:
                pass
            else:
                fraction=float(line.split()[4])
                acceptor=line.split()[0]
                donor1=line.split()[1]
                donor2=line.split()[2]
                res1_chain, res1_name, orig_res1_num, atom1=utilities.parse_hbond(acceptor, reverse_dict, accept=True)
                res2_chain, res2_name, orig_res2_num, atom2=utilities.parse_hbond(donor1, reverse_dict)
                if fraction*100 < 1:
                    break
                hbond='%s.%s%s@%s-%s.%s%s@%s' % (res1_chain, res1_name,orig_res1_num,atom1,res2_chain, res2_name, orig_res2_num,atom2)
                percent=fraction*100
                print hbond, 'persists for %0.2f %% of simulation' % percent
        fhandle.close()
    return

def clustering(cwd, outname, ref, trjfile, cluster, d=2.0):
    ref=os.path.abspath(ref)
    ref_basename=os.path.basename(ref)
    trjfile=os.path.abspath(trjfile)
    #make sure have absolute path since working in analysis folder now
    make_analysis_folder(cwd, 'clustering')
    #check if already clustered with this outname
    if not os.path.exists('%s-d%s-cluster-centroid-pdbs/' % (outname, d)):
        os.mkdir('%s-d%s-cluster-centroid-pdbs/' % (outname,d))
    else:
        print "Existing clustering results for %s-d%s-cluster-centroid-pdbs/ so will not overwite" % (outname, d)
        sys.exit()
    
    # test that the residues you cluster on are the right ones
    new_ambmask=utilities.check_cluster_pdbfile(ref, cluster, outname)
    external_file_io.write_clustering_ptraj_file(ref, trjfile, new_ambmask, d, outname)
    command='%s/bin/cpptraj %s ptraj-cluster.in' %  (os.environ['AMBERHOME'], ref)
    output, err=run_linux_process(command)
    numpy.savetxt('ptraj-cluster.out', output.split('\n'), fmt='%s')
    if 'rror' in err:
        numpy.savetxt('ptraj-cluster.err', err.split('\n'), fmt='%s')
        print "ERROR IN CPPTRAJ RMSD CALCULATION"
        print "CHECK ptraj-cluster.err"
        print err
        sys.exit()
    fhandle=open('%s-cluster-d%s-summary_out' % (outname, d))
    count=0
    num=0
    total_fraction=0
    for line in fhandle.readlines():
        if '#' in line:
            pass
        else:
            num+=1
            fraction=float(line.split()[2])
            total_fraction+=fraction
            if total_fraction > 0.90:
                print "%s clusters represent 90 percent of simulation" % num
                break
    for x in range(0, (num)):
        shutil.move('%s-cluster-d%s-centroid.c%s.pdb' % (outname, d, x), '%s-d%s-cluster-centroid-pdbs/' % (outname, d))
    extra_files=glob.glob('%s-cluster-d%s-centroid*pdb' % (outname, d))
    for file in extra_files:
        os.remove(file)
    return

def hbonds(cwd, outname, ref, trjfile, selection):
    if len(selection) != 2:
        print "HBOND REQUIRES 2 INPUT GROUPS"
        sys.exit()
    ref=os.path.abspath(ref)
    ref_basename=os.path.basename(ref)
    trjfile=os.path.abspath(trjfile)
    make_analysis_folder(cwd, 'hbonds')
    #make sure have absolute path since working in analysis folder now
    residue_mapper=utilities.map_residues(ref)
    new_ambmask1=utilities.check_cluster_pdbfile(ref, [selection[0],], outname)
    new_ambmask2=utilities.check_cluster_pdbfile(ref, [selection[1],], outname)
    external_file_io.write_hbond_ptraj(trjfile, new_ambmask1, new_ambmask2, outname)
    command='%s/bin/cpptraj %s ptraj-hbonds.in' %  (os.environ['AMBERHOME'], ref)
    output, err=run_linux_process(command)
    numpy.savetxt('ptraj-hbonds.out', output.split('\n'), fmt='%s')
    if 'rror' in err:
        numpy.savetxt('ptraj-hbonds.err', err.split('\n'), fmt='%s')
        print "ERROR IN CPPTRAJ RMSD CALCULATION"
        print "CHECK ptraj-cluster.err"
        print err
        sys.exit()
    files=glob.glob('avghb-%s-*.dat' % outname)
    external_file_io.write_all_hbonds_to_pml(files, outname, residue_mapper,  ref)
    print "wrote %s_hbonds.pml" % outname
    return


def selection_checker(cwd, reffile, selection):
    ref=os.path.abspath(reffile)
    ref_basename=os.path.basename(ref)
    make_analysis_folder(cwd, 'selection')
    residue_mapper=utilities.map_residues(ref)
    new_ambmask=utilities.check_cluster_pdbfile(ref, selection, outname='test')
    return



def main(args):
    if args.debug==True:
        import pdb
        pdb.set_trace()
    try: 
        os.environ['AMBERHOME']
    except KeyError:
        print "AMBERHOME environment variable is not set"
        print "On AWS this is /home/mlawrenz/amber14/"
        sys.exit()
    cwd=os.getcwd()
    if args.reffile is None:
        print "SUPPLY A REFERENCE PDBFILE"
        sys.exit()
    if not os.path.exists(args.reffile):
        print "SUPPLY A REFERENCE PDBFILE"
        sys.exit()
    if args.analysis=='selection_check':
        selection_checker(cwd, args.reffile, args.selection)
        sys.exit()
    if args.trjfile is None:
        print "SUPPLY A TRJFILE"
        sys.exit()
    if not os.path.exists(args.trjfile):
        print "SUPPLY A TRJFILE"
        sys.exit()
    if args.analysis=='rmsd_calc':
        print "Running %s analysis" % args.rmsdtype
        rmsd_and_rmsf(cwd, args.reffile, args.trjfile, args.rmsdtype)
    if args.analysis=='clustering':
        print "Running clustering analysis" 
        clustering(cwd, args.outname, args.reffile, args.trjfile, args.selection, args.distance)
    if args.analysis=='hbonds':
        hbonds(cwd, args.outname, args.reffile, args.trjfile, args.selection)
    if args.analysis=='solvent_calc':    
        sovlent_calc(cwd, args.outname, args.reffile, args.trjfile, args.selection, args.occupancy)
    print "done with analysis"
    


def parse_commandline():
    parser = argparse.ArgumentParser(description='''
Run chosen analysis on MD DCD trajectory file. Requires a reference PDBfile and
DCD trajectory passed in. Choose analysis workflow: rmsd or clustering. See additional options for each
analysis choice by running: 
$SCHRODINGER/run /home/mlawrenz/AnalyzeMD/AnalyzeMD.py clustering -h
Examples:
$SCHRODINGER/run ~/AnalyzeMD/AnalyzeMD.py -r reference.pdb -t trj.dcd  rmsd_calc rmsd-all
$SCHRODINGER/run ~/AnalyzeMD/AnalyzeMD.py -r reference.pdb -t trj.dcd clustering --selection '50-55,131.B' ':48,49.F' -o inter-chain -d 1.0
$SCHRODINGER/run ~/AnalyzeMD/AnalyzeMD.py  -r reference.pdb -t trj.dcd hbonds -o 4mdk-holo-redo --group1 1-163 --group2 164-240^C

''')


    #parser.add_argument("analysis", choices=["rmsd", "clustering", "hbond"])
    parser.add_argument('-r', '--reffile', dest='reffile',
                      help='reference structure PDB file, for RMSD and visualization. MUST MATCH TRAJECTORY', required=True)
    parser.add_argument('-t', '--trjfile', dest='trjfile', help='DCD trj for trajectory analysis')
    parser.add_argument('--debug', dest='debug', action="store_true")
    subparsers= parser.add_subparsers(help='analysis suboption choices', dest='analysis')
    mask_parser=argparse.ArgumentParser(add_help=False)
    outname_parser=argparse.ArgumentParser(add_help=False)
    outname_parser.add_argument('--outname', dest='outname', default='protein',
                      help='name for hbonds output, please do not use dashes or weird characters')
    mask_parser.add_argument('--selection', nargs='*', dest='selection', help='''Residues to use for analysis, with separated by
commas and dashes and chains specified at the end with a period. If you do not
provide a chain then we assign it to chain A.  i.e. 45-55.A 68,128,155.B period.
Can include multiple selections from multiple chains.
''')
    x_parser=subparsers.add_parser("selection_check", parents=[mask_parser], help='''
pass in the mask and get the corresponding AMBER numbers''')
    r_parser=subparsers.add_parser("rmsd_calc", help='''
Computes RMSD (distance from reference structure) and RMSF (distance from average
structure) for all heavy atoms or just backbone atoms. All output is written as
described in the README section of rmsd-analysis/ directory that is created.
The chosen datatype is mapped to the bfactors of reference file for visualization.
bfactor-${choice}-${referencefile}.pdb
''')
    r_parser.add_argument('rmsdtype',choices=['rmsd-all', 'rmsd-bb', 'rmsf-all', 'rmsf-bb'], help='type of rmsd or rmsf calculation, using all heavy atoms or just backbone CA atoms')
    c_parser=subparsers.add_parser("clustering", parents=[mask_parser], help='''
Uses hierarchical clustering and average linkage for RMSD of selected residues
passed in with selection, and the minimum distance between clusters is
specified by the distance argument.''')
    c_parser.add_argument('-d', '--distance', dest='distance', default=2.0,
                      help='clustering RMSD cutoff to define clusters. recommend 2 for most small changes.')
    c_parser.add_argument('-o', '--outname', dest='outname', default='set',
                      help='name for clustering output, will precede a name NAME-cluster-*')
    h_parser=subparsers.add_parser("hbonds", parents=[mask_parser], help='''
Compute hydrogen bonds between two groups of residues specified with --selection
(can be all to all with same selection. Output is pml file with "strong" hbonds
defined as  >= 50.0 and red, "medium" >= 10 and < 50.0 and pink, and
"weak" as >=1 and < 10 and yellow.''')
    h_parser.add_argument('-o', '--outname', dest='outname', default='set',
                      help='name for hbonds output, please do not use dashes or weird characters')
    c_parser=subparsers.add_parser("solvent_calc", parents=[mask_parser, outname_parser], help='''
Analyze water occupancy over trajectory and compute solvent hbonds to a provided
mask. Will output a PDB file with water density at a %% of the max density.
Default this is 0.8. Also will report significant solvent hbonds.''')
    c_parser.add_argument('-o', '--occupancy', dest='occupancy', default=0.8,
                      help='Occupancy >=X%% of the max density name used to output PDB file with these waters.')
    args = parser.parse_args()
    return args

#run the function main if namespace is main
if __name__ == "__main__":
    args = parse_commandline()
    main(args)
