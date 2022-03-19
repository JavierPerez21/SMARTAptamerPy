import os, sys
from pyfastaq import sequences
import argparse
import re
import random

def random_nucleotide(match):
    return random.sample(["A", "C", "G", "T"], k=1)[0]

def cleanfastq(indir, outdir=None):
    # Rel range refers to the relevant indices of the sequences for the study
    # Modify this line to fit your needs
    rel_range = [24, 54]

    if outdir is None:
        outdir = indir[:-4]

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    i = 0
    for infile in os.listdir(indir):
        newfile = outdir + "/" + str(infile)
        if infile.endswith(".fastq"):
            with open(indir + "/" + infile, "r") as inf:
                lines = inf.readlines()
                with open(newfile, "w") as outf:
                    for i in range(0, len(lines), 4):
                        data = lines[i:i+4]
                        if data[1][rel_range[0]:rel_range[1]].count("N") > 2:
                            continue
                        else:
                            data[1] = re.sub("N", random_nucleotide, data[1])
                            outf.write("".join(data))
                outf.close()
            inf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--molecule', type=str, default="target",
                        help='Molecule name to be used as target')
    parser.add_argument('-r', '--r', type=str, default="R1",
                        help='R1 for forward sequence R2 for comlementary')
    args = parser.parse_args()
    molecule = args.molecule
    r = args.r
    indir = str(os.getcwd()) + f"/data/{molecule}_fastq_{r}"
    cleanfastq(indir)