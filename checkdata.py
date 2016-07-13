import numpy
import pylab

data=numpy.loadtxt('time-hb-Ub-ligand-accept.dat', skiprows=1, dtype=int)
file=open('time-hb-Ub-ligand-accept.dat')
header=file.readline()
header=header.split()
file.close()
for x in range(0,10):
    print header[x]
    print len(numpy.where(data[:,x]!=0)[0])/float(len(data[:,x]))

import itertools
import mutinf
combos=itertools.combinations(range(1,10), 2)
mis=dict()
corrs=dict()
for combo in combos:
    x=combo[0]
    y=combo[1]
    mi=mutinf.mutual_information_2d(data[:,x], data[:,y])
    mis[combo]=mi
    r=numpy.corrcoef(data[:,x], data[:,y])
    corrs[combo]=r[0][1]

import operator
sorted_corrs=sorted(corrs.iteritems(), key=operator.itemgetter(1))
sorted_mis=sorted(mis.iteritems(), key=operator.itemgetter(1))
for term in reversed(sorted_corrs):
    x=term[0][0]
    y=term[0][1]
    corr=term[1]
    if numpy.abs(corr) < 0.2:
        #print header[x], header[y], 'R=%s' % corr, 'MI=%s' % mis[(x,y)]
        pylab.figure()
        pylab.scatter(range(0, len(data[:,x])), data[:,x], color='r', label=header[x])
        pylab.title(header[x])
        pylab.figure()
        pylab.scatter(range(0, len(data[:,x])), data[:,y], color='k', label=header[x])
        pylab.title(header[y])
        print header[x], header[y], 'R=%s' % corr
        pylab.show()

