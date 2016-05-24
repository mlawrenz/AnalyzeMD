import glob
import argparse
import numpy
import os
import sys

def percent_score(percent):
    if percent >= 50.0:
        color='red'
    if percent >= 10 and percent < 50.0:
        color='pink' 
    if percent < 10:
        color='yellow'
    return color

def parse_all_hbonds(files, prefix, residues, ref):
    format_data=dict()
    all_hbonds=[]
    pymol_handle=open('%s_pymol_hbond_label.py' % prefix, 'w')
    pymol_handle.write('load %s\n' % ref)
    pymol_handle.write('sel protein, polymer\n')
    pymol_handle.write('sel lhr, resi 357-367\n')
    pymol_handle.write('sel ring, resi 380-427\n')
    pymol_handle.write('show cartoon, protein\n')
    pymol_handle.write('hide lines, protein\n')
    ohandle=open('all_hbonds_%s.dat' % prefix, 'w')
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
            res1_num=int(acceptor.split('@')[0].split('_')[1])   
            res1_name=acceptor.split('@')[0].split('_')[0]
            orig_res1_num=resid[int(res1_num)]
            atom1=acceptor.split('@')[1].split('-')[0]
            res2_num=int(donor.split('@')[0].split('_')[1])   
            res2_name=donor.split('@')[0].split('_')[0]
            orig_res2_num=resid[int(res2_num)]
            atom2=donor.split('@')[1].split('-')[0]
            hbond='%s%s@%s-%s%s@%s' % (res1_name,orig_res1_num,atom1,res2_name, orig_res2_num,atom2)
            print hbond, percent
            ohandle.write('%s\t%s\n' % (hbond, percent))
            res1_chain='A'
            res2_chain='A'
            pymol_accept='%s/%s/%s' % (res1_chain, orig_res1_num, atom1)
            pymol_donor='%s/%s/%s' % (res2_chain, orig_res2_num, atom2)
            pymol_handle.write('show sticks, chain %s and resi %s\n' % (res1_chain, orig_res1_num))
            pymol_handle.write('show sticks, chain %s and resi %s\n' % (res2_chain, orig_res2_num))
            hcolor=percent_score(percent)
            if hcolor=='red':
                pymol_handle.write('dist %s_strong_hbonds, %s, %s\n' % (prefix, pymol_donor, pymol_accept))
                pymol_handle.write('color %s, %s_strong_hbonds\n' % (hcolor, prefix))
            if hcolor=='pink':
                pymol_handle.write('dist %s_medium_hbonds, %s, %s\n' % (prefix, pymol_donor, pymol_accept))
                pymol_handle.write('color %s, %s_medium_hbonds\n' % (hcolor, prefix))
            if hcolor=='yellow':
                pymol_handle.write('dist %s_weak_hbonds, %s, %s\n' % (prefix, pymol_donor, pymol_accept))
                pymol_handle.write('color %s, %s_weak_hbonds\n' % (hcolor, prefix))
    pymol_handle.close()
    ohandle.close()
    return 


