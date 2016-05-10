#!/usr/bin/env python
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


def map_file_to_resnum(datafile, pdb_data):
    """
    Parses input from a two column data file, creating a dictionary of
    column[0]:column[1] key/value pairs for resid:data correspondence.  
    Goes through pdb line by line.  If residue is in dictionary data_dict,
    the b-factor is replaced in output_file by the dictionary value.  If the
    residue is not in the dictionary, the b-factor is given value 0.0.
    Returns void.
    """
    # Read in file
    col1=[]
    col2=[]
    tmp=open(datafile)
    for line in tmp:
        if '#' in line:
            continue
        col1.append(float(line.split()[0]))
        col2.append(float(line.split()[1]))
    #import numpy
    #col1=numpy.loadtxt(datafile, usecols=(0,))
    #col2=numpy.loadtxt(datafile, usecols=(1,))

    data_dict = dict()
    for (resid, d) in zip(col1, col2):
        data_dict[int(resid)]=d


    out = []
    for line in pdb_data:
        if line[0:6] == "ATOM  ":
            resnum = int(line[23:26].strip())
            if resnum in data_dict.keys():
                out.append("%s%6.2F%s" % (line[:60],data_dict[resnum],
                                          line[66:]))
            else:
                out.append("%s%6s%s" % (line[:60],"NA",line[66:]))
        
        elif line[0:6] == "HETATM":
            out.append("%s%6s%s" % (line[:60],"NA",line[66:]))
        
        else:
            out.append(line)
    
    return out

def map_file_by_index(datafile, pdb_data):
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
    print "in map function"
    for line in pdb_data:
        if line[0:6] == "ATOM  ":
            resnum = int(line[23:26].strip())
            if prev_resnum==0:
                prev_resnum=resnum
            if resnum!=prev_resnum:
                prev_resnum=resnum
                index+=1
            out.append("%s%6.2F%s" % (line[:60],col2[index],
                                          line[66:]))
        elif line[0:6] == "HETATM":
            print "not changing for HETATM entry"
            out.append("%s%6s%s" % (line[:60],"NA",line[66:]))
        
        else:
            out.append(line)
    
    return out


def main(pdbfile, datafile, mapresid=False):
#    """
#    Call if program called from command line.
#    """
    f = open(pdbfile,'r')
    pdb_data = f.readlines()
    f.close()

    if mapresid==True:
        #print "Mapping data by resid from datafile to pdbfile"
        #print "Make sure atoms correpond!"
        print "THIS ISN'T WORKING R N"
        #out=map_file_to_resnum(datafile, pdb_data)
    else:
        print "running file map"
        out=map_file_by_index(datafile, pdb_data)
        print "Warning: Mapping data by residue index from datafile to pdbfile"
        print "Warning: Make sure atoms correpond!"
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

