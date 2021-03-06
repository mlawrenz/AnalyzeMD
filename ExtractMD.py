from schrodinger.trajectory.desmondsimulation import create_simulation
import numpy
import external_file_io
import shlex
import shutil
import subprocess
from subprocess import PIPE
import optparse
from schrodinger.structure import write_ct, write_ct_pdb, StructureReader
from optparse import OptionParser
from schrodinger.structutils.analyze import evaluate_asl
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.application.desmond.util import get_indices, parse_slice
import os
from trajectory_extract_utilities import *


"""
Extract MD

"""


def _parse_frames(frames_str):
    slice_list=[]
    for s in frames_str.split(','):
        s = s.strip()
        if s:
            slice_list.append(parse_slice(s, inclusive=True))
    return slice_list


def run_linux_process(command):
    #command=command.replace("//", "/")
    #command_list=shlex.split(command)
    p=subprocess.Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    p.wait()
    output, err=p.communicate()
    return output, err

def write_vmd_script(cmsfile, trjfile, step=1, ligand=False):
    basename=os.path.basename(trjfile).split('_trj')[0]
    ohandle=open('wrap.tcl', 'w')
    ohandle.write(''' 
mol new {0} type mae first 0 last -1 step {2} waitfor all
mol addfile {1}/clickme.dtr type dtr first 0 last -1 step 5 waitfor all
package require pbctools
set cell [ pbc get -now ]
pbc set $cell'''.format(cmsfile, trjfile, step))
    if ligand!=False: 
        ohandle.write('''
set solute [ atomselect top "protein or resname {0}" ]')
pbc wrap -centersel "resname {0}" -center com -compound residue -all'''.format(ligand))
    else:
        ohandle.write('''
set solute [ atomselect top "protein" ]'
pbc wrap -centersel "protein" -center com -compound residue -all''') 
    ohandle.write('''
set n [molinfo top get numframes]
for {{set t 0}} {{$t < $n}} {{incr t}} {{
set all [atomselect top "all" frame $t]
set ref [atomselect top "backbone" frame 0]
set bb [atomselect top "backbone" frame $t]
set trans_mat [measure fit $bb $ref]
$all move $trans_mat
}}
animate write dcd aln-wrap-{0}.dcd sel $all waitfor all 
animate goto start
animate write pdb aln-wrap-{0}.pdb sel $all beg 0 end 1
exit'''.format(basename))
    return

def extract_traj_with_water(cwd, cmsfile, trjfile, asl_expr, step=1):
    cmsfile=os.path.abspath(cmsfile)
    trjfile=os.path.abspath(trjfile)
    basename=os.path.basename(trjfile).split('_trj')[0]
    if asl_expr != 'protein':
        st = StructureReader(cmsfile).next()
        asl_searcher = AslLigandSearcher()
        ligands = asl_searcher.search(st)
        if len(ligands) != 1:
            print "Only one ligand must meet definition"
            sys.exit()
        else:
            res=ligands[0].pdbres
        write_vmd_script(cmsfile, trjfile, step, res)
    else:
        write_vmd_script(cmsfile, trjfile, step)
    command='%s/vmd -dispdev text -e wrap.tcl' % (os.environ['VMD_DIR'])
    # AGAIN, POpen is not working!
    #os.system(command)
    output, err=run_linux_process(command)
    if 'rror' in output:
        numpy.savetxt('vmd.err', output.split('\n'), fmt='%s')
        print "ERROR IN VMD"
        print "CHECK vmd.err"
        print output.split('\n')
        sys.exit()
    if not os.path.exists('aln-wrap-%s.dcd' % basename):
        print "DID NOT WRITE DCD FILE"
        sys.exit()
    if not os.path.exists('%s/analysis' % cwd):
        os.mkdir('%s/analysis' % cwd)
    shutil.move('aln-wrap-%s.dcd' % basename,  '%s/analysis' % cwd)
    shutil.move('aln-wrap-%s.pdb' % basename,  '%s/analysis' % cwd)
    print "trajectories are in analysis/ directory"
    return

def extract_frames_from_trajectory(options):
    #create a temporary directory to run the model
    cmsfile=os.path.abspath(options.cmsfile)
    trjfile=os.path.abspath(options.trjfile)
    tmp_folder = "/tmp/extractmd_%s/" % (os.getpid())
    if not os.path.exists(tmp_folder) :
        os.mkdir(tmp_folder)
    os.chdir(tmp_folder)

    asl_expr = options.asl
    basename=os.path.basename(trjfile).split('_trj')[0]
    csim = create_simulation(cmsfile, trjfile)

    # hardcode basename
    basename = 'solute-%s' % basename
    frames='::5'
    total_frame = csim.total_frame
    slice_list = _parse_frames(frames)
    frames = get_indices(slice_list, total_frame)
    # only extract PDBs
    extract_func = extract_pdb
    for i in frames:
        filename = "%s/%s_%05d.%s"%(tmp_folder, basename, i, 'pdb')
        extract_func(csim, i, filename, asl = asl_expr )
    print "wrote %s pdbfiles" % i
    print "trajectories are in analysis/ directory"
    return tmp_folder, basename


def convert_desmond_to_dcd(cwd, tmp_folder, basename):
    program='%s/catdcd' % os.environ['CATDCD_DIR']
    if not os.path.exists('%s/analysis' % cwd):
        os.mkdir('%s/analysis' % cwd)
    command='{0} -o {1}/analysis/{2}.dcd -otype dcd -s {3}/{2}_00000.pdb -pdb `ls -v {3}/{2}*pdb`'.format(program, cwd, basename, tmp_folder)
    # RUNNING THIS WAY BC CATCHING OUTPUT DID NOT WORK FOR SOME REASON
    os.system(command)
    command='cp {0}/{1}_00000.pdb {2}/analysis/{1}_reference.pdb'.format(tmp_folder, basename, cwd)
    output, err=run_linux_process(command)
    if 'rror' in err:
        numpy.savetxt('copy.err', err.split('\n'), fmt='%s')
        print "ERROR IN COPY"
        print "CHECK copy.err"
        print err
        sys.exit()
    return

def main(options):
    cwd=os.getcwd()
    try:
        os.environ['AMBERHOME']
    except KeyError:
        print "AMBERHOME environment variable is not set"
        print "On AWS this is /home/mlawrenz/amber14/"
        sys.exit()
    try:
        os.environ['CATDCD_DIR']    
    except KeyError:
        print "CATDCD_DIR environment variable is not set"
        print "On AWS this is /home/mlawrenz/linuxamd64/bin/catdcd4.0/"
        sys.exit()
    if options.debug==True:
        import pdb
        pdb.set_trace()
    if options.solvent==True:
        print "Extracting wrapped DCD trajectory with solvent"
        print "Warning: may still have jumps in protein, use only for water analysis"
        try:
            os.environ['VMD_DIR']
        except KeyError:
            print "VMD_DIR environment variable is not set"
            print "On AWS this is /home/mlawrenz/VMD1.9.2/bin/"
            sys.exit()
        extract_traj_with_water(cwd, options.cmsfile, options.trjfile, options.asl, step=1)
    else:
        print "This may take up to tens of minutes"
        tmp_folder, basename=extract_frames_from_trajectory(options)
        print "Converting to solvent-free DCD trajectory"
        convert_desmond_to_dcd(cwd, tmp_folder, basename)
        os.chdir('%s/analysis' % cwd)
        external_file_io.write_ptraj_align_file('%s_reference.pdb' % basename ,'%s.dcd' % basename)
        command='%s/bin/cpptraj %s ptraj-align.in' %  (os.environ['AMBERHOME'], '%s_reference.pdb' % basename)
        output, err=run_linux_process(command)
        if 'rror' in err:
            numpy.savetxt('ptraj-align.err', err.split('\n'), fmt='%s')
            print "ERROR IN CPPTRAJ ALIGN CALCULATION"
            print "CHECK ptraj-align.err"
            print err
            sys.exit()
        os.remove('%s.dcd' % basename)
        print "Cleaning up"
        shutil.rmtree(tmp_folder)
    return
    


def parse_commandline():
    parser = optparse.OptionParser(description='''Converts a single Schrodinger Desmond
trajectory to a single DCD trajectory with a corresponding PDB reference file. This
exact reference file needs to be used for analysis, with the residue numbers and
chain information. Use AnalyzeMD.py to align this PDB file to a reference crystal PDB, and combine
trajectories from multiple runs into a single trajectory.
Required environment variables (as installed on AWS):                                                                                 
export AMBERHOME=/home/mlawrenz/amber14                                                              
export CATDCD_DIR=/home/mlawrenz/linuxamd64/bin/catdcd4.0/                                          
export VMD_DIR=/home/mlawrenz/VMD1.9.2/bin/
''')

    parser.add_option('-c', '--cmsfile', dest='cmsfile',
                      help='desmond cmsfile for trajectory to analyze')
    parser.add_option('-t', '--trjfile', dest='trjfile',
                      help='desmond trj/ directory for trajectory analysis')
    parser.add_option("-a", "--asl",
                  action="store", dest="asl", default='protein',
                  help='asl to specify what part of structure to be extracted. Default is protein. Can use "protein or ligand"')
    parser.add_option('--solvent', action="store_true",  help='Extract solvent trajs. These may have a few remaining jumps in the protein. Use only for solvent analysis (which still may have some problems).')
    parser.add_option('--debug', dest='debug', action="store_true")
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(options)


