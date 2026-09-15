#Usage: python cpmdidx.py -c data/cpmd.inp -t data/topol.top resid is 107 and name is HD2

import argparse
import mimicpy

parser = argparse.ArgumentParser()
parser.add_argument("sele", nargs='+', help="Selection query", type=str)
parser.add_argument("-c", "--cpmd", help="CPMD filename", default='cpmd.inp', type=str)
parser.add_argument("-t", "--top", help="MPT filename", default='topol.mpt', type=str)
parser.set_defaults(feature=False)
args = parser.parse_args()

sele = ' '.join(args.sele)
cpmd_file = args.cpmd
mpt_file = args.top

cpmd = mimicpy.CpmdScript.from_file(cpmd_file)
ids = {int(i.split()[1]):int(i.split()[3]) for i in cpmd.MIMIC.OVERLAPS.splitlines()[1:]}
topol = mimicpy.Mpt.from_file(mpt_file)
selection = topol.select(sele)
print("The atoms are:")
for i in selection.index:
    print(ids[i])
