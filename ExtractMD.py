from schrodinger.trajectory.desmondsimulation import create_simulation
import shutil
import subprocess
from subprocess import PIPE
import optparse
import schrodinger.trajectory.analysis as analysis
from schrodinger.structure import write_ct, write_ct_pdb
from optparse import OptionParser
from schrodinger.structutils.analyze import evaluate_asl
from schrodinger.application.desmond.util import get_indices, parse_slice
import os
from trajectory_extract_utilities import *


"""
Analyze MD
===========

** To be used on Schrodinger Desmond results  **

Program requires as input (see options with -h flag):
*
*

Program requires AMBERTools and will produce DCD files
This file can be loaded into VMD:
* you can run the command: vmd -pdb reference.pdb -dcd file.dcd

"""

def _parse_frames(frames_str):
    slice_list=[]
    for s in frames_str.split(','):
        s = s.strip()
        if s:
            slice_list.append(parse_slice(s, inclusive=True))
    return slice_list


def run_linux_process(command):
    p=subprocess.Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    p.wait()
    output, err=p.communicate()
    return output, err

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
    basename = 'extract-%s' % basename
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
    return tmp_folder, basename


def convert_desmond_to_dcd(cwd, tmp_folder, basename):
    program='/home/mlawrenz/LINUXAMD64/bin/catdcd4.0/catdcd'
    if not os.path.exists('%s/analysis' % cwd):
        os.mkdir('%s/analysis' % cwd)
    command='{0} -o {1}/analysis/{2}.dcd -otype dcd -s {3}/{2}_00000.pdb -pdb `ls -v {3}/{2}*pdb`'.format(program, cwd, basename, tmp_folder)
    os.system(command)
    #output, err=run_linux_process(command)
    #print output.split('\n')
    #if 'rror' in err:
    #    numpy.savetxt('catdcd.err', err.split('\n'), fmt='%s')
    #    print "ERROR IN CATDCD"
    #    print "CHECK catdcd.err"
    #    print err
    #    sys.exit()
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
    tmp_folder, basename=extract_frames_from_trajectory(options)
    print "Converting to DCD"
    convert_desmond_to_dcd(cwd, tmp_folder, basename)
    print "Cleaning up"
    shutil.rmtree(tmp_folder)
    return
    


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-c', '--cmsfile', dest='cmsfile',
                      help='cmsfile for trajectory to analyze')
    parser.add_option('-t', '--trjfile', dest='trjfile',
                      help='trj directory for trajectory analysis')
    parser.add_option("-a", "--asl",
                  action="store", dest="asl", default='protein',
                  help="asl to specify what part of structure to be extracted")
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(options)


