import glob
from schrodinger.structutils.analyze import evaluate_asl
from schrodinger.structure import write_ct, write_ct_pdb, StructureReader
import external_file_io
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
Utility functions
"""


def run_linux_process(command):
    p=subprocess.Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    p.wait()
    output, err=p.communicate()
    return output, err



def get_residues_from_radius(radius, ref):
    asl = 'fillres ((within %i ligand ) and not (res.ptype "T3P"))' % (int(float(radius)))
    st=StructureReader(ref).next()
    atoms=evaluate_asl(st, asl)
    if len(atoms)==0:
        print "NO LIGAND IN REFERENCE FILE"
        sys.exit()
    pocket=st.extract(atoms)
    resnames=[i.pdbres for i in pocket.residue]
    resnums=[i.resnum for i in pocket.residue]
    chains=[i.chain for i in pocket.residue]
    data=dict()
    for chain in set(chains):
        for (n,c) in zip(resnums, chains):
            if c not in data.keys():
                data[c]=[]
            if str(n) not in data[c]:
                data[c].append(str(n))
    selection=[]
    for chain in data.keys():
        selection.append(','.join(data[chain]) + '.%s' % chain)
    if len(selection)==1:
        selection=[selection,]
    return selection



def parse_hbond(amber_mask, reverse_dict, accept=False):
    if 'Solvent' in amber_mask:
        chain='water'
        orig_num=''
        name=''
        atom='H'
        if accept==True:
            atom='O'
    else:
        num=amber_mask.split('@')[0].split('_')[1]
        name=amber_mask.split('@')[0].split('_')[0]
        chain=reverse_dict[num][0]
        orig_num=reverse_dict[num][1]
        atom=amber_mask.split('@')[1].split('-')[0]
    return chain, name, orig_num, atom


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
        if line.split()[0]=='ATOM' or line.split()[0]=='HETATM':
            if line.split()[3] in exclude:
                break
            if 'pseu' in line:
                break
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
    external_file_io.write_strip_ptraj_file(ref, new_ambmask, outname)
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


def make_keypaths(nested):
    for key, value in nested.iteritems():
        if isinstance(value, collections.Mapping):
            for subkey, subvalue in make_keypaths(value):
                yield [key] + subkey, subvalue
        else:
            yield [key], value


