import glob
import optparse
from optparse import OptionParser
import external_file_io
import utilities
import collections
import shutil
import numpy
import sys
import pdb_bfactor
import os
import subprocess
from subprocess import PIPE
import argparse
import textwrap as _textwrap

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

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width,
initial_indent=indent, subsequent_indent=indent) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

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
        bb=pdb_bfactor.check_pdb_file(pdb_data)
        maxval, out=pdb_bfactor.map_file_by_index(datafile, pdb_data, bb)
    else:
        maxval, out=pdb_bfactor.map_file_by_index(datafile, pdb_data)
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
    for line in output.split('\n'):
        if 'Coordinate processing' in line:
            total_frames=int(line.split()[5])
            print "processed %s frames" % total_frames
    if selection!=None:
        solvent_amber_data=utilities.process_solvent_output(outname, reverse_dict)
        file='%s_bridgeout.dat' % outname
        external_file_io.write_bridged_residues_to_pml(file, solvent_amber_data, total_frames, outname, reverse_dict, ref)
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
        print "Your selection is:", selection
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


def traj_combine(cwd, outname, file_list):
    program='%s/catdcd' % os.environ['CATDCD_DIR']
    command='{0} -o {1}/combine-{2}.dcd -otype dcd -dcd `cat {3}`'.format(program, cwd, outname, file_list)
    output, err=run_linux_process(command)
    erroroutput=False
    if 'rror' in err:
        erroroutput=err
    if 'rror' in output:
        erroroutput=output
    if erroroutput==True: 
        numpy.savetxt('catdcd.err', erroroutput.split('\n'), fmt='%s')
        print "ERROR IN CATDCD CALCULATION"
        print "CHECK catdcd.err"
        sys.exit()
    if not os.path.exists('{0}/combine-{1}.dcd'.format(cwd, outname)):
        print "ERROR IN CATDCD CALCULATION"
    else:
        print "combined dcd trajetcory is {0}/combine-{1}.dcd".format(cwd, outname)
    return

def pdb_align(ref, inputfile, asl='backbone'):
    verdict=utilities.check_selection_boolean(ref, 'T3P')
    if verdict==True:
        print "FOR NOW DO NOT ALIGN SOLVENT STRUCTURES THIS WAY"
        sys.exit()
    verdict=utilities.check_selection_boolean(inputfile, 'T3P')
    if verdict==True:
        print "FOR NOW DO NOT ALIGN SOLVENT STRUCTURES THIS WAY"
        sys.exit()
    basename_file=os.path.basename(inputfile)
    program='%s/utilities/structalign' % os.environ['SCHRODINGER']
    command='{0} -asl "{1}" {2} {3}'.format(program, asl, ref, inputfile)
    output, err=run_linux_process(command)
    erroroutput=False
    if 'rror' in err:
        erroroutput=err
    if 'rror' in output:
        erroroutput=output
    if erroroutput==True: 
        numpy.savetxt('align.err', erroroutput.split('\n'), fmt='%s')
        print "ERROR IN ALIGN CALCULATION"
        print "CHECK align.err"
        sys.exit()
    if not os.path.exists('rot-%s' % basename_file):
        print "ERROR IN ALIGN CALCULATION"
    else:
        shutil.move('rot-%s' % basename_file, 'aligned-%s' % basename_file)
        print '%s is aligned to %s, named aligned-%s' % (inputfile, ref, basename_file)
        for line in output.split('\n'):
            if 'RMSD' in line:
                value=float(line.split()[1])
                if value > 1.0:
                    print "WARNING, RMSD FIT IS HIGH at %0.2f" % value
                    print "BE SURE TO CONFIRM YOUR ALIGNMENT IS CORRECT"
                    print "TRY PASSING IN --radius around ligand, if present"
    return

def main(args):
    cwd=os.getcwd()
    if args.debug==True:
        import pdb
        pdb.set_trace()
    # run jobs that don't need other checks
    if args.analysis=='traj_combine':
        try:
            os.environ['CATDCD_DIR']
            traj_combine(cwd, args.outname, args.file_list)
            sys.exit()
        except KeyError:
            print "CATDCD_DIR environment variable is not set"
            print "On AWS this is /home/mlawrenz/linuxamd64/bin/catdcd4.0/"
            sys.exit()
    try:
        os.environ['SCHRODINGER']
    except KeyError:
        print "SCHRODINGER environment variable is not set"
        sys.exit()
    if args.reffile is None:
        print "SUPPLY A REFERENCE PDBFILE"
        sys.exit()
    if not os.path.exists(args.reffile):
        print "SUPPLY A REFERENCE PDBFILE"
        sys.exit()
    if args.analysis=='pdb_align':
        if args.radius!=None:
            asl='fillres within 10 (ligand)'
            pdb_align(args.reffile, args.inputfile, asl)
        else:
            pdb_align(args.reffile, args.inputfile)
        sys.exit()
    # run jobs that need AMBER
    try: 
        os.environ['AMBERHOME']
    except KeyError:
        print "AMBERHOME environment variable is not set"
        print "On AWS this is /home/mlawrenz/amber14/"
        sys.exit()
    # check for traj file
    if args.trjfile is None:
        print "SUPPLY A TRJFILE"
        sys.exit()
    if not os.path.exists(args.trjfile):
        print "SUPPLY A TRJFILE"
        sys.exit()
    if args.analysis=='rmsd_calc':
        print "Running %s analysis" % args.rmsdtype
        rmsd_and_rmsf(cwd, args.reffile, args.trjfile, args.rmsdtype)
        sys.exit()
    # all analysis below requires these
    if args.radius is None and args.selection is None:
        print "Need to pass in a selection or radius around ligand"
        sys.exit()
    if args.radius is not None:
        selection=utilities.get_residues_from_radius(args.radius, args.reffile)
    else:
        selection=args.selection
    if args.analysis=='selection_check':
        selection_checker(cwd, args.reffile, selection)
        sys.exit()
    if args.analysis=='clustering':
        print "Running clustering analysis" 
        clustering(cwd, args.outname, args.reffile, args.trjfile, selection, args.distance)
    if args.analysis=='hbonds':
        hbonds(cwd, args.outname, args.reffile, args.trjfile, selection)
    if args.analysis=='solvent_calc':    
        sovlent_calc(cwd, args.outname, args.reffile, args.trjfile, selection, args.occupancy)
    print "done with analysis"
    


def parse_commandline():
    parser = argparse.ArgumentParser(description='''
Requires a reference PDBfile (all residues numbers determined from this) and DCD
trajectory. *MAKE SURE YOU ALIGN THIS FILE TO YOUR PROJECT NEEDS*
 Choose analysis from
positional arguments (NOTE THAT POSITIONS FOR ARGS MATTER). See additional options for each analysis choice by running: 
|n
$SCHRODINGER/run /home/mlawrenz/AnalyzeMD/AnalyzeMD.py clustering -h
|n
Examples:
|n
$SCHRODINGER/run ~/AnalyzeMD/AnalyzeMD.py -r reference.pdb -t trj.dcd  rmsd_calc rmsd-all
|n
$SCHRODINGER/run ~/AnalyzeMD/AnalyzeMD.py -r reference.pdb -t trj.dcd clustering --selection '50-55,131.B' '48,49.F' -o inter-chain -d 1.0
$SCHRODINGER/run ~/AnalyzeMD/AnalyzeMD.py -r reference.pdb -t trj.dcd clustering --radius 4 -o protein-ligand -d 1.0
|n
$SCHRODINGER/run ~/AnalyzeMD/AnalyzeMD.py  -r reference.pdb -t trj.dcd hbonds -o inter-chain --selection '*.B' '*.C'
$SCHRODINGER/run ~/AnalyzeMD/AnalyzeMD.py  -r reference.pdb -t trj.dcd hbonds -o protein-ligand --radius 5
|n
$SCHRODINGER/run ~/AnalyzeMD/AnalyzeMD.py  -r waters-reference.pdb -t waters-trj.dcd solvent_calc  -o 'ligand' --selection 'DRG.L'
''', usage='%(prog)s [OPTIONS]', formatter_class=MultilineFormatter)


    #parser.add_argument("analysis", choices=["rmsd", "clustering", "hbond"])
    parser.add_argument('-r', '--reffile', dest='reffile',
                      help='reference structure PDB file, for RMSD and visualization. MUST MATCH TRAJECTORY', required=True)
    parser.add_argument('-t', '--trjfile', dest='trjfile', help='DCD trj for trajectory analysis')
    parser.add_argument('--debug', dest='debug', action="store_true")
    subparsers= parser.add_subparsers(help='analysis suboption choices', dest='analysis')
    mask_parser=argparse.ArgumentParser(add_help=False)
    outname_parser=argparse.ArgumentParser(add_help=False)
    outname_parser.add_argument('-o', '--outname', dest='outname', default='protein', required=True, help='name for output, please do not use dashes or weird characters')
    mask_parser.add_argument('--selection', nargs='*', dest='selection', help='''Residues to use for analysis, with separated by
commas and dashes and chains specified at the end with a period. If you do not
provide a chain then we assign it to chain A.  i.e. 45-55.A 68,128,155.B period.
Can include multiple selections from multiple chains.
''')
    radius_parser=argparse.ArgumentParser(add_help=False)
    radius_parser.add_argument('--radius',  dest='radius', help='''Radius around ligand to define selection of residues for analysis''')
    a_parser=subparsers.add_parser("pdb_align", parents=[radius_parser], help='Align 2 structure files')
    a_parser.add_argument('-f', '--inputfile', dest='inputfile',
                      help='file to align to your reference passed in with -r')
    d_parser=subparsers.add_parser("traj_combine", parents=[outname_parser], help='Combine DCD files into one')
    d_parser.add_argument('-f', '--file-list', dest='file_list',
                      help='file with DCD trajectories to be combined')
    x_parser=subparsers.add_parser("selection_check", parents=[mask_parser, radius_parser], help='''
pass in the mask and get the corresponding AMBER numbers''')
    r_parser=subparsers.add_parser("rmsd_calc", help='''
Computes RMSD (distance from reference structure) and RMSF (distance from average
structure) for all heavy atoms or just backbone atoms. All output is written as
described in the README section of rmsd-analysis/ directory that is created.
The chosen datatype is mapped to the bfactors of reference file for visualization.
bfactor-${choice}-${referencefile}.pdb
''')
    r_parser.add_argument('rmsdtype',choices=['rmsd-all', 'rmsd-bb', 'rmsf-all', 'rmsf-bb'], help='type of rmsd or rmsf calculation, using all heavy atoms or just backbone CA atoms')
    c_parser=subparsers.add_parser("clustering", parents=[mask_parser, radius_parser, outname_parser], help=''' Uses hierarchical clustering and average linkage for RMSD of selected residues
passed in with selection, and the minimum distance between clusters is
specified by the distance argument.''')
    c_parser.add_argument('-d', '--distance', dest='distance', default=2.0,
                      help='clustering RMSD cutoff to define clusters. recommend 2 for most small changes.')
    h_parser=subparsers.add_parser("hbonds", parents=[mask_parser, radius_parser, outname_parser], help=''' Compute hydrogen bonds between two groups of residues specified with --selection
(can be all to all with same selection. Output is pml file with "strong" hbonds
defined as  >= 50.0 and red, "medium" >= 10 and < 50.0 and pink, and
"weak" as >=1 and < 10 and yellow.''')
    c_parser=subparsers.add_parser("solvent_calc", parents=[mask_parser, radius_parser, outname_parser], help=''' Analyze water occupancy over trajectory and compute solvent hbonds to a provided
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
