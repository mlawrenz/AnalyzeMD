import numpy
import pylab
import glob
import itertools
import operator
import utilities


def correlations(dir, residue_mapper, plot=False):
    reverse_dict = {value: keypath for keypath, value in utilities.make_keypaths(residue_mapper)}
    files=glob.glob('%s/time*dat' % dir)
    n=0
    for file in files:
        tmp_data=numpy.loadtxt(file, skiprows=1, dtype=int)
        f=open(file)
        tmp_header=f.readline()
        f.close()
        if n==0:
            data=tmp_data
            amber_header=tmp_header.split()
            n+=1
        else:
            import pdb
            pdb.set_trace()
            for column in range(1, tmp_data.shape[1]):
                new=tmp_data[:,column].reshape(tmp_data[:,column].shape[0], 1)
                data=numpy.hstack((data, new))
                amber_header.append(tmp_header.split()[column])
    header=[]
    for h in amber_header:
        if 'Frame' in h:
        header.append(h)
            pass
        else:
            first=h.split('-')[0]
            second=h.split('-')[0]
            res1_chain, res1_name, orig_res1_num, atom1=utilities.parse_hbond(first, reverse_dict)
            res2_chain, res2_name, orig_res2_num, atom2=utilities.parse_hbond(second, reverse_dict)
            hbond='%s.%s%s@%s-%s.%s%s@%s' % (res1_chain, res1_name,orig_res1_num,atom1,res2_chain, res2_name, orig_res2_num,atom2)
            header.append(hbond)

    combos=itertools.combinations(range(1,10), 2)
    corrs=dict()
    for combo in combos:
        x=combo[0]
        y=combo[1]
        r=numpy.corrcoef(data[:,x], data[:,y])
        corrs[combo]=r[0][1]
    sorted_corrs=sorted(corrs.iteritems(), key=operator.itemgetter(1))
    for term in reversed(sorted_corrs):
        x=term[0][0]
        y=term[0][1]
        corr=term[1]
        if numpy.abs(corr) > 0.2:
            pylab.figure()
            pylab.scatter(range(0, len(data[:,x])), data[:,x], color='r', label=header[x])
            pylab.title(header[x])
            pylab.figure()
            pylab.scatter(range(0, len(data[:,x])), data[:,y], color='k', label=header[x])
            pylab.title(header[y])
            print header[x], header[y], 'R=%s' % corr
            if plot==True:
                pylab.show()
    return

