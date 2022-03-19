import os, sys
import argparse
from pyfastaq import sequences
import re
import random
def random_nucleotide(match):
    return random.sample(["A", "C", "G", "T"], k=1)[0]

def fastaq2txt(molecule='theophylline', indir=None, r='r1'):
    # Rel range refers to the relevant indices of the sequences for the study
    # r1 and r2 refers to the sequence of study: r1 for the normal sequences and r2 for the complementary
    # Modify these lines to fit your needs
    if r == 'r1':
        rel_range = [24, 54]
    elif r == 'r2':
        rel_range = [18, 48]


    if indir is None:
        outdir = str(os.getcwd()) + "/data/{}_txt_{}/".format(molecule, r)
    else:
        outdir = indir.replace("fastq", "txt")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for infile in os.listdir(indir):
        if infile.endswith(".fastq"):
            seqs = []
            seq_reader = sequences.file_reader(indir + "/" + infile)
            for sequence in seq_reader:
                seq = sequence.seq[rel_range[0]: rel_range[1]]
                if seq.count("N") < 3:
                    seq = re.sub("N", random_nucleotide, seq)
                    seqs.append(seq)
            outfile = infile.replace(".fastq", ".txt")
            with open(outdir + "/" + outfile, 'w') as fout:
                for seq in seqs:
                    fout.write(seq + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--molecule', type=str, default="target",
                        help='Molecule name to be used as target')
    parser.add_argument('-r', '--r', type=str, default="R1",
                        help='R1 for forward sequence R2 for comlementary')
    args = parser.parse_args()
    molecule = args.molecule
    r = args.r
    indir = str(os.getcwd()) + "/data/{}_fastq_{}".format(molecule, r)
    fastaq2txt(molecule, indir=indir, r=r)
