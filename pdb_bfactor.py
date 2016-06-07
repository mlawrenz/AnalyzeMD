#!/usr/bin/env python
import numpy
import optparse
import sys

__description__ = \
"""
pdb_bfactor.py:

Alters the bfactor column of a pdb file. 
Default assigns the data by index to the ATOM entries of a PDBfile.
Optional approach maps the resid in the 2 column data file to the resid in the
PDBfile.
to read into pymol use:
    spectrum b, blue_white_red, minimum=20, maximum=50
    as cartoon
    cartoon putty
"""


def check_pdb_file(pdb_data):
    prev_resnum=0
    index=0
    bb=[]
    compare=dict()
    for line in pdb_data:
        if 'ATOM' in line or 'HETATM' in line:
            resnum = int(line[23:26].strip())
            if line.split()[2]=='CA':
                bb.append(resnum)
    return bb


def map_file_by_index(datafile, pdb_data, bb=False):
#    """
#    adds first data to first residue, without resid mapping
#    col2=
    col1=[]
    col2=[]
    tmp=open(datafile)
    for line in tmp:
        if '#' in line:
            continue
        col1.append(float(line.split()[0]))
        col2.append(float(line.split()[1]))

    out = []
    prev_resnum=0
    index=0
    #print "in map function"
    check=[]
    num_residues=1
    filestart=True
    for line in pdb_data:
        if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
            resnum = int(line[23:26].strip())
            if filestart==True:
                prev_resnum=resnum
                check.append('%s%s' % (resnum, line.split()[3]))
                filestart=False
            if resnum!=prev_resnum:
                prev_resnum=resnum
                check.append('%s%s' % (resnum, line.split()[3]))
                if bb is not False:
                    if not resnum in bb:
                        print "excluding entry ", line.split()[3]
                    else:
                        num_residues+=1
                        index+=1
                        pass
                else:
                   num_residues+=1
                   index+=1
            if bb is not False:
                if not resnum in bb:
                    out.append("%s%s%s" % (line[:60],'  NA', line[66:]))
                else:
                    out.append("%s%6.2F%s" % (line[:60],col2[index], line[66:]))
            else:
                out.append("%s%6.2F%s" % (line[:60],col2[index], line[66:]))
        else:
            out.append(line)
    
    if len(col1) != num_residues:
        print "MISMATCH IN NUMBER OF ATOMS"
        print "PDBFILE %s AND DATAFILE %s" % (num_residues, len(col1))
        sys.exit()
    else:
        print "Wrote data for %s residues to bfactor-pdbfile" % num_residues
    # pass this as estimate for the pymol color spectrum
    maxval=numpy.percentile(col2, 90)
    return maxval, out


def main(pdbfile, datafile, bb_flag=False):
#    """
#    Call if program called from command line.
#    """
    f = open(pdbfile,'r')
    pdb_data = f.readlines()
    f.close()
    if bb_flag==True:
        bb=check_pdb_file(pdb_data)
        maxval, out=map_file_by_index(datafile, pdb_data, bb)
    else:
        print "Mapping data to PDB file"
        maxval, out=map_file_by_index(datafile, pdb_data)
    outfilename='bfactor-%s' % pdbfile
    outfile=open(outfilename, 'w')
    for line in out:
        outfile.write(line)
    outfile.close()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-p', '--pdbfile', type='string', dest='pdbfile',
                      help='pdbfile')
    parser.add_option('-d', '--datafile', type='string', dest='datafile',
                      help='2 column data with B factor column info')
    parser.add_option('-r', '--mapresid', dest='mapresid', action="store_true")
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    if options.mapresid==True:
        main(options.pdbfile, options.datafile, mapresid=True)
    else:
        main(options.pdbfile, options.datafile)

