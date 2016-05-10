import sys
import os
from schrodinger.trajectory.desmondsimulation import create_simulation
import schrodinger.trajectory.analysis as analysis
from schrodinger.structure import write_ct, write_ct_pdb
from optparse import OptionParser
from schrodinger.structutils.analyze import evaluate_asl
from schrodinger.application.desmond.util import get_indices, parse_slice

def _parse_frames(frames_str):
    slice_list=[]
    for s in frames_str.split(','):
        s = s.strip()
        if s:
            slice_list.append(parse_slice(s, inclusive=True))
    return slice_list

def extract_ct(csim, frameno, filename, asl=None):
    st = csim.getFrameStructure(frameno)
    if asl:
        selected_indices = evaluate_asl(st, asl)
        st = st.extract(selected_indices, copy_props=True)

    # delete obsolete properties
    properties = ['s_chorus_trajectory_file', 's_m_original_cms_file', 's_ffio_ct_type']
    for prop in properties:
        if prop in st.property:
            del st.property[prop]

    return st

def extract_mae(csim, frameno, filename, asl=None):
    st = extract_ct(csim, frameno, filename, asl=asl)
    write_ct(st, filename)

def extract_pdb(csim, frameno, filename, asl=None):
    st = extract_ct(csim, frameno, filename, asl=asl)
    write_ct_pdb(st, filename)

def extract_cms(csim, frameno, filename, asl=None):
    if asl:
        raise NotImplementedError
    frame = csim.getFrame(frameno)
    csim.cst.writeCms(filename, frame)

if __name__ == "__main__":

    class MyParser(OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    usage = """$SCHRODINGER/run -FROM desmond trajectory_extract_frame.py <cmsfile> <trjfile>
A script to extract selected frames of trajectory into a series of output structure files."""
    epilog = '''
----------------Examples-------------------------------------------------------

Example1: Extract all frames and write each frame into maestro file
  $SCHRODINGER/run -FROM desmond trajectory_extract_frame.py \\
                    input-out.cms input_trj

Example2: Extract just the protein from your system and fifth frame.
          Use basename 'prot'
          NOTE: The frame selection format is "from:to:skip"
                empty field suggests from=0; to=last; skip=1
  $SCHRODINGER/run -FROM desmond trajectory_extract_frame.py \\
                    input-out.cms input_trj -a 'protein' -f '::5' -b 'prot'

Example3: Extract frames 10, 20-40, and frame 71. Output in pdb format.
  $SCHRODINGER/run -FROM desmond trajectory_extract_frame.py \\
                    input-out.cms input_trj -f '10,20:40,71' -o pdb
'''

    parser = MyParser(usage=usage,epilog=epilog)

    parser.add_option("-a", "--asl",
                  action="store", dest="asl", default=None,
                  help="asl to specify what part of structure to be extracted")

    parser.add_option("-l", "--asl_file",
                  action="store", dest="asl_file", default=None,
                  help="asl from a file to specify what part of structure to be extracted")

    parser.add_option("-b", "--basename",
                  action="store", dest="basename", default=None,
                  help="basename of the output structure")

    parser.add_option("-o", "--output",
                  action="store", dest="output", default='mae',
                  choices = ['mae', 'cms', 'pdb'],
                  help="format of the output structure")

    parser.add_option("-f", "--frames",
                  action="store", dest="frames", default=None,
                  help="inclusive range expressions to specify frames. "       +
                       "For instance, ':4:2, 5' selects frame 0, 2, 4 and 5. " +
                       "Note that negative integer number can also be used in "+
                       "the range expression.""")

    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 2:
        parser.print_help()
        sys.exit(1)

    if options.asl and options.asl_file:
        print "--asl and --asl_file are exclusive."
        sys.exit(1)

    asl_expr = None
    if options.asl:
        asl_expr = options.asl
    elif options.asl_file:
        try:
            asl_expr = ''.join(open(options.asl_file).readlines())
            asl_expr = asl_expr.strip()
        except Exception, e:
            print "Fail to parse ASL from %s: %s"%(options.asl_file, e)
            sys.exit(1)

    cmsfile = args[0]
    trjfile = args[1]
    csim = create_simulation(cmsfile, trjfile)

    basename = options.basename
    if not basename:
        base_filename = os.path.basename(cmsfile)
        t = base_filename.split('.')
        if len(t) >= 2:
            basename = '.'.join(t[:-1])
        else:
            # cmsfile does not contain suffix
            basename = base_filename

    total_frame = csim.total_frame
    if options.frames == None:
        frames = range(total_frame)
    else:
        slice_list = _parse_frames(options.frames)
        frames = get_indices(slice_list, total_frame)
    if options.output == 'mae':
        extract_func = extract_mae
    elif options.output == 'cms':
        extract_func = extract_cms
        if options.asl:
            print "-a options is not supported when output structure format is .cms."
            sys.exit(1)
    elif options.output == 'pdb':
        extract_func = extract_pdb

    for i in frames:
        filename = "%s_%05d.%s"%(basename, i, options.output)
        extract_func(csim, i, filename, asl = asl_expr )

