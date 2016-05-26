import glob
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
cluster {2} mass epsilon {3} out {4}-cluster-d{3}_out averagelinkage gracecolor summary {4}-cluster-d{3}-summary_out info {4}-cluster-d{3}-Cluster_info repout {4}-cluster-d{3}-centroid repfmt pdb
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

def parse_hbond(amber_mask, reverse_dict):
    num=amber_mask.split('@')[0].split('_')[1]
    name=amber_mask.split('@')[0].split('_')[0]
    # mapper values are amber
    if num not in reverse_dict.keys():
        # probably doing a solvent hbond
        chain='water'
        orig_num=num
    else:
        chain=reverse_dict[num][0]
        orig_num=reverse_dict[num][1]
    atom=amber_mask.split('@')[1].split('-')[0]
    return chain, name, orig_num, atom


def write_solvent_ptraj_file(ref, trjfile, mask, outname, occupancy):
#hbond H1 donormask :1-64.O= acceptormask :1-64.O= nointramol solventdonor :WAT solventacceptor :WAT.O o ut numhb.dat avgout avghb.dat solvout avgsolvent.dat bridgeout bridge.dat 
#radial rad.dat 0.1 10.0 :T3P@O :346
#radial rad2.dat 0.1 10.0 :T3P@H1 :346 
#radial rad3.dat 0.1 10.0 :T3P@H2 :346 
    ohandle=open('ptraj-water.in', 'w')
    ohandle.write('''
trajin {0}
reference {1}
# Create average of solute to view with grid.
center {2}
rms reference @C,CA,O,N
grid {3}_water.dx 50 0.5 50 0.5 50 0.5 :T3P@O pdb {3}_water{4}.pdb max {4}
hbond H1 donormask {2}@N*,O* acceptormask :T3P&@N*,O* nointramol out numhb-ligand-donor.dat  avgout avghb-ligand-donor.dat series 
hbond H2 acceptormask {2}&@N*,O* donormask :T3P&@N*,O* nointramol out numhb-ligand-accept.dat avgout avghb-ligand-accept.dat series 
'''.format(trjfile, ref, mask, outname, occupancy))
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

def make_analysis_folder(cwd, name):
    if not os.path.exists('%s/%s-analysis' % (cwd, name)):
        os.mkdir('%s/%s-analysis' % (cwd, name))
    os.chdir('%s/%s-analysis' % (cwd, name))
    return

def catch_output_and_errors(output, error, name):
    ohandle=open('%s.err' % name, 'w')
    for l in output:
        ohandle.write('%s\n' % l)
    ohandle.close()
    
    numpy.savetxt('ptraj-rmsd.err', err, fmt='%s')

def parse_ambmask(residue_mapper, selection):
    total_residues_list=[]
    new_residues_list=[]
    for mask in selection:
        chain=mask.split('.')[-1]
        if '*' in mask:
            if chain=='*':
                for chain in residue_mapper.keys():
                    for res in sorted(residue_mapper[chain].keys()):
                        total_residues_list.append('%s.%s' % (chain, str(res)))
                        new_residues_list.append(residue_mapper[chain][res])
            else:
                #chain is specified
                for res in sorted(residue_mapper[chain].keys()):
                    total_residues_list.append('%s.%s' % (chain, str(res)))
                    new_residues_list.append(residue_mapper[chain][res])
        else:
            residues_list=mask.split('.')[0].split(',')
            for x in residues_list:
                x=x.strip(':')
                if '-' in x:
                    start=int(x.split('-')[0])
                    end=int(x.split('-')[1])
                    for i in range(start, end+1):
                        total_residues_list.append('%s.%s' % (chain, str(i)))
                        new_residues_list.append(residue_mapper[chain][str(i)])
                else:
                    total_residues_list.append('%s.%s' % (chain, str(x)))
                    new_residues_list.append(residue_mapper[chain][str(x)])
    return total_residues_list, new_residues_list

def map_residues(ref):
    exclude=['T3P', 'NA ', 'Cl ']
    residue_mapper=dict()
    pdbfile=open(ref)
    amber_resnum=1
    prev_resnum=0
    for line in pdbfile.readlines():
        if line.split()[3] in exclude:
            break
        if 'pseu' in line:
            break
        if line.split()[0]=='ATOM' or line.split()[0]=='HETATM':
            resnum = str(line[23:26].strip())
            chain = str(line[20:22].strip())
            if chain not in residue_mapper.keys():
                residue_mapper[chain]=dict()
            if prev_resnum==0:
                prev_resnum=resnum
                residue_mapper[chain][resnum]=str(amber_resnum)
                amber_resnum+=1
            if resnum!=prev_resnum:
                prev_resnum=resnum
                residue_mapper[chain][resnum]=str(amber_resnum)
                amber_resnum+=1
        else:
            pass
    return residue_mapper


        
def check_cluster_pdbfile(ref, selection, outname):
    residue_mapper=map_residues(ref)
    total_residue_list, new_residue_list=parse_ambmask(residue_mapper, selection)
    new_ambmask=','.join([str(i) for i in new_residue_list])
    new_ambmask=':%s' % new_ambmask
    write_strip_ptraj_file(ref, new_ambmask, outname)
    #write_strip_ptraj_file(ref, ambmask, outname)
    command='%s/bin/cpptraj %s ptraj-strip.in' %  (os.environ['AMBERHOME'], ref)
    output, err=run_linux_process(command)
    pdbfile=open('%s-cluster-check.pdb' % outname)
    amber_residues=[]
    prev_resnum=0
    for line in pdbfile.readlines():
        if 'ANISOU' in line:
            print "PDB reference contains ANISOU lines, please remove these first"
            sys.exit()
        if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
            resnum = int(line[23:26].strip())
            if prev_resnum==0:
                amber_residues.append(''.join([line.split()[3], line.split()[5]]))
                prev_resnum=resnum
            elif resnum!=prev_resnum:
                prev_resnum=resnum
                amber_residues.append(''.join([line.split()[3], line.split()[5]]))
            else:
                pass
    pdbfile.close()
    print "You asked for %s" % selection
    #print "We converted to %s" % new_ambmask
    print "You are getting:"
    print amber_residues
    return new_ambmask

def percent_score(percent):
    if percent >= 50.0:
        color='red'
    if percent >= 10 and percent < 50.0:
        color='pink' 
    if percent < 10:
        color='yellow'
    return color


def keypaths(nested):
    for key, value in nested.iteritems():
        if isinstance(value, collections.Mapping):
            for subkey, subvalue in keypaths(value):
                yield [key] + subkey, subvalue
        else:
            yield [key], value



def parse_all_hbonds_to_pml(files, outname, residue_mapper, ref):
    format_data=dict()
    all_hbonds=[]
    pymol_handle=open('%s_hbonds.pml' % outname, 'w')
    pymol_handle.write('load %s\n' % ref)
    pymol_handle.write('sel protein, polymer\n')
    pymol_handle.write('show cartoon, protein\n')
    pymol_handle.write('hide lines, protein\n')
    ohandle=open('%s_all_hbonds.dat' % outname, 'w')
    # need a damn reverse dict for this
    reverse_dict = {value: keypath for keypath, value in keypaths(residue_mapper)}
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
            res1_chain, res1_name, orig_res1_num, atom1=parse_hbond(acceptor, reverse_dict)
            res2_chain, res2_name, orig_res2_num, atom2=parse_hbond(donor, reverse_dict)
            hbond='%s.%s%s@%s-%s.%s%s@%s' % (res1_chain, res1_name,orig_res1_num,atom1,res2_chain, res2_name, orig_res2_num,atom2)
            print hbond, percent
            ohandle.write('%s\t%s\n' % (hbond, percent))
            pymol_accept='%s/%s/%s' % (res1_chain, orig_res1_num, atom1)
            pymol_donor='%s/%s/%s' % (res2_chain, orig_res2_num, atom2)
            pymol_handle.write('show sticks, chain %s and resi %s\n' % (res1_chain, orig_res1_num))
            pymol_handle.write('show sticks, chain %s and resi %s\n' % (res2_chain, orig_res2_num))
            hcolor=percent_score(percent)
            if hcolor=='red':
                pymol_handle.write('dist %s_strong_hbonds, %s, %s\n' % (outname, pymol_donor, pymol_accept))
                pymol_handle.write('color %s, %s_strong_hbonds\n' % (hcolor, outname))
            if hcolor=='pink':
                pymol_handle.write('dist %s_medium_hbonds, %s, %s\n' % (outname, pymol_donor, pymol_accept))
                pymol_handle.write('color %s, %s_medium_hbonds\n' % (hcolor, outname))
            if hcolor=='yellow':
                pymol_handle.write('dist %s_weak_hbonds, %s, %s\n' % (outname, pymol_donor, pymol_accept))
                pymol_handle.write('color %s, %s_weak_hbonds\n' % (hcolor, outname))
    pymol_handle.write('hide labels, ligand_weak_hbonds\n')
    pymol_handle.write('hide labels, ligand_medium_hbonds\n')
    pymol_handle.write('hide labels, ligand_strong_hbonds\n')
    pymol_handle.close()
    ohandle.close()
    return 


def rmsd_and_rmsf(cwd, ref, trjfile, datatype):
    ref=os.path.abspath(ref)
    ref_basename=os.path.basename(ref)
    trjfile=os.path.abspath(trjfile)
    #make sure have absolute path since working in analysis folder now
    make_analysis_folder(cwd, 'rmsd')
    write_rmsd_ptraj_file(ref, trjfile)
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
    write_load_bfactor_pml(datatype, ref_basename, maxval)
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

def sovlent_calc(cwd, outname, ref, trjfile, selection, occupancy=0.8):
    ref=os.path.abspath(ref)
    ref_basename=os.path.basename(ref)
    trjfile=os.path.abspath(trjfile)
    residue_mapper=map_residues(ref)
    reverse_dict = {value: keypath for keypath, value in keypaths(residue_mapper)}
    #make sure have absolute path since working in analysis folder now
    make_analysis_folder(cwd, 'solvent')
    # test that the residues you cluster on are the right ones
    new_ambmask=check_cluster_pdbfile(ref, selection, outname)
    write_solvent_ptraj_file(ref, trjfile, new_ambmask, outname, occupancy)
    command='%s/bin/cpptraj %s ptraj-water.in' %  (os.environ['AMBERHOME'], ref)
    output, err=run_linux_process(command)
    if 'rror' in err:
        numpy.savetxt('ptraj-cluster.err', err.split('\n'), fmt='%s')
        print "ERROR IN CPPTRAJ RMSD CALCULATION"
        print "CHECK ptraj-cluster.err"
        print err
        sys.exit()
    for file in glob.glob('avghb-*.dat'):
        fhandle=open(file)
    count=0
    num=0
    total_fraction=0
    for line in fhandle.readlines():
        if '#' in line:
            pass
        else:
            fraction=float(line.split()[4])
            acceptor=line.split()[0]
            donor1=line.split()[1]
            donor2=line.split()[2]
            res1_chain, res1_name, orig_res1_num, atom1=parse_hbond(acceptor, reverse_dict)
            res2_chain, res2_name, orig_res2_num, atom2=parse_hbond(donor1, reverse_dict)
            if fraction*100 < 1:
                break
            hbond='%s.%s%s@%s-%s.%s%s@%s' % (res1_chain, res1_name,orig_res1_num,atom1,res2_chain, res2_name, orig_res2_num,atom2)
            percent=fraction*100
            print hbond, 'persists for %0.2f %% of simulation' % percent
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
    new_ambmask=check_cluster_pdbfile(ref, cluster, outname)
    write_clustering_ptraj_file(ref, trjfile, new_ambmask, d, outname)
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
    residue_mapper=map_residues(ref)
    new_ambmask1=check_cluster_pdbfile(ref, [selection[0],], outname)
    new_ambmask2=check_cluster_pdbfile(ref, [selection[1],], outname)
    write_hbond_ptraj(trjfile, new_ambmask1, new_ambmask2, outname)
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
    parse_all_hbonds_to_pml(files, outname, residue_mapper,  ref)
    print "wrote %s_hbonds.pml" % outname
    return


def selection_checker(cwd, reffile, selection):
    ref=os.path.abspath(reffile)
    ref_basename=os.path.basename(ref)
    make_analysis_folder(cwd, 'selection')
    residue_mapper=map_residues(ref)
    new_ambmask=check_cluster_pdbfile(ref, selection, outname='test')
    return



def main(args):
    if args.debug==True:
        import pdb
        pdb.set_trace()
    try: 
        os.environ['AMBERHOME']
    except KeyError:
        print "AMBERHOME environment variable is not set"
        print "On AWS this is /home/mlawrenz/VMD1.9.2/bin/"
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




