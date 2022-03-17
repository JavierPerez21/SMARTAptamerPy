import argparse
import os
import multiprocessing
import subprocess
import re
import numpy as np
import pandas as pd
import random
import datetime


class MDARNAFolder():
    def __init__(self, k=6, threads=multiprocessing.cpu_count(), temp=37,
                 output_dir=str(os.getcwd()) + "/outputs",
                 ):
        self.k = k
        self.threads = threads
        self.temp = temp
        self.output_dir = output_dir

    def run(self):
        self.get_g4_structure()
        self.process_g4_structure()
        self.calculate_mda()
        subprocess.call([f"rm -rf {self.output_dir + '/mfold_seeds'}"], shell=True)   # Why is this not working?
        subprocess.call([f"rm -rf {self.output_dir + '/structure'}"], shell=True)     # Why is this not working?

    def get_g4_structure(self):
        print("Getting g4 structure...")
        structure_dict = self.output_dir + '/structure'
        if not os.path.exists(structure_dict):
            os.makedirs(structure_dict)
        seedsdir = self.output_dir + '/mfold_seeds'
        if not os.path.exists(seedsdir):
            os.makedirs(seedsdir)
        subprocess.call([f"split -n l/{self.threads} -d --additional-suffix=mfold {self.output_dir}/represent.fasta {seedsdir}/"], shell=True)
        tseeds = [seed for seed in os.listdir(seedsdir) if seed.endswith('mfold')]
        #tseeds = tseeds.split()

        pm = multiprocessing.Pool(self.threads)

        for file in tseeds:
            G4 = {}
            ids = []
            with open(seedsdir + "/" + file, 'r') as inf1:
                for line1 in inf1:
                    info = line1.strip().split('\t')
                    ids.append(info[0])
                    id = info[0]
                    with open(self.output_dir + '/structure/' + id + '.txt', 'w') as outf:
                        outf.write(info[1] + '\n')
                    outf.close()
                    bool = subprocess.check_output([f"./qgrs-cpp/qgrs -i {self.output_dir}/structure/{id}.txt -o {self.output_dir}/structure/{id}.seq.txt"], shell=True)
                    if bool.split() == ["No", "QGRS", "found"]:
                        G4[id] = 0
                    else:
                        with open(self.output_dir + '/structure/' + id + '.seq.txt', 'r') as inf2:
                            for line in inf2:
                                if line.strip() == '':
                                    continue
                                elif line.strip().startswith('ID'):
                                    continue
                                elif line.strip().startswith('-'):
                                    continue
                                else:
                                    result = line.strip().split()
                                    G4[id] = result[6]
                        inf2.close()
            inf1.close()
            with open(seedsdir + '/' + file + '-g4.txt', 'w') as outf:
                for id in ids:
                    outf.write(id + '\t')
                    if id in G4.keys():
                        outf.write(str(G4[id]) + '\n')
                    else:
                        outf.write('0\n')
            outf.close()
        pm.close()
        pm.join()

        subprocess.call([f"cat {seedsdir}/*-g4.txt > {self.output_dir}/g4.txt"], shell=True)

    def process_g4_structure(self):
        print("Processing g4 structure...")
        duplicates = {}
        s2i = {}
        with open(self.output_dir + '/represent.fasta', 'r') as inf:
            with open(self.output_dir + '/structure/seq.txt', 'w') as outf:
                for line in inf:
                    info = line.strip().split('\t')
                    if info[1] not in s2i:
                        s2i[info[1]] = info[0]
                    else:
                        if info[1] not in duplicates:
                            duplicates[info[1]] = [info[0]]
                        else:
                            duplicates[info[1]].append(info[0])
                    outf.write(info[1] + '\n')
        inf.close()
        outf.close()

        rand = random.randint(0, 5000)
        cmd = f"RNAfold --jobs={self.threads} --infile={self.output_dir}/structure/seq.txt --outfile=structure-log_{str(rand)}.txt --noPS -T {self.temp} --noconv"
        subprocess.call([cmd], shell=True)
        subprocess.call([f"mv {str(os.getcwd())}/structure-log_{str(rand)}.txt {self.output_dir}/structure-log.txt"], shell=True)
        print(cmd)

        seq = ''
        dg = {}
        outf = open(self.output_dir + '/structure.txt', 'w')
        with open(self.output_dir + '/structure-log.txt', 'r') as inf:
                for line in inf:
                    if re.search('[ATCG]+', line):
                        seq = line.strip()
                    elif re.search('.*\((.*?)\)', line):
                        dg[s2i[seq]] = line.strip().split('(')[-1].split(')')[0]
        inf.close()

        for seq in duplicates:
            for key in duplicates[seq]:
                dg[key] = dg[s2i[seq]]

        g4 = {}
        with open(self.output_dir + '/g4.txt', 'r') as inf:
            for line in inf:
                info = line.strip().split('\t')
                g4[info[0]] = info[1]
        inf.close()

        print(dg.keys())
        print(g4.keys())
        for key in g4.keys():
            outf.write(key + '\t' + dg[key] + '\t' + g4[key] + '\n')
        outf.close()

    def calculate_rscore(self, seq, kmer_score):
        lkmer = self.k
        nu = 0
        score = 0
        scores = []
        while nu <= (len(seq) - lkmer):
            kmer = seq[nu:nu + lkmer]
            if kmer in kmer_score.keys():
                score += kmer_score[kmer]
            nu += 1
        score = score / 10
        return score

    def calculate_mda(self):
        """
        Calculated MDA score
        :return:
        """
        print("Calculating MDA score...")
        family_size = {}
        with open(self.output_dir + '/aptamer_clusters', 'r') as inf:
            for line in inf:
                info = line.strip().split('\t')
                family_size[info[0]] = int(info[1])
        inf.close()

        # recaculate score based on conserved score
        kmer_score = {}
        with open(self.output_dir + '/score_kmers.txt', 'r') as inf:
            for line in inf:
                info = line.strip().split('\t')
                kmer_score[info[0]] = float(info[1])
        inf.close()

        # calculate rscore
        rscore = {}
        sequence = {}
        with open(self.output_dir + '/represent.fasta', 'r') as inf:
            for line in inf:
                info = line.strip().split('\t')
                rscore[info[0]] = self.calculate_rscore(info[1], kmer_score)
                sequence[info[0]] = info[1]
        inf.close()

        scores = sorted(rscore.values())
        max_score = scores[-1]
        min_score = scores[0]
        for key in rscore.keys():
            rscore[key] = (rscore[key] - min_score) * 10 / (max_score - min_score)

        # calculate structure score
        structure_score = {}
        dg = {}
        g4 = {}
        with open(self.output_dir + '/structure.txt', 'r') as inf:
            for line in inf:
                if re.search('dG', line):
                    pass
                else:
                    info = line.strip().split('\t')
                    dg[info[0]] = -float(info[1]) / 2
                    g4[info[0]] = float(info[2]) / 9.5
        inf.close()
        scores = sorted(dg.values())
        max_score = scores[-1]
        min_score = scores[0]
        for key in dg.keys():
            # dg[key] = (dg[key] - min_score) * 10 / (max_score - min_score)   # why are we doing this here, feels like an error done on purpose
            if dg[key] > g4[key]:
                structure_score[key] = dg[key]
            else:
                structure_score[key] = g4[key]

        # calculate score for family size
        family_score = {}
        sizes = list(family_size.values())
        mean = int(np.mean(sizes)) + 1
        sizes = sorted(sizes)
        max = sizes[-1]
        print(max, mean)
        for key in family_size.keys():
            nu = 10
            if max > mean * nu * 10:
                while nu >= 0:
                    if family_size[key] > mean * nu * 10:
                        family_score[key] = nu
                        break
                    nu -= 1
            else:
                while nu >= 0:
                    if family_size[key] > mean * nu:
                        family_score[key] = nu
                        break
                    nu -= 1

        # calculate final score
        final_score = {}
        presentation = {}
        for key in rscore.keys():
            fkey = key.split('_seq_')[0]
            scores = [family_score[fkey], rscore[key], structure_score[key]]
            scores = sorted(scores)
            final_score[key] = scores[0] / 2 + scores[1] / 2
            if scores[2] == family_score[fkey]:
                presentation[key] = 'KS-SS'
            elif scores[2] == rscore[key]:
                presentation[key] = 'FS-SS'
            elif scores[2] == structure_score[key]:
                presentation[key] = 'FS-KS'

        data = []
        nu = 1
        with open(self.output_dir + '/result_aptamer.txt', 'w') as outf:
            outf.write('Rank\tfamily_sequence_id\tFS:family_size|family_size_score\tKS:kmer_score\tSS:dG|qgrs-G4score|structure_score\tscore\tscore_based\tsequence\n')
            for key in sorted(final_score, key=final_score.get, reverse=True):
                fkey = key.split('_seq_')[0]
                outf.write(f"{nu}\t{fkey}\t{family_size[fkey]}\|{family_score[fkey]}\t{rscore[key]}\t{dg[key]}\|{g4[key]}\|{structure_score[key]}\t{final_score[key]}\t{presentation[key]}\t{sequence[key]}\n")
                data.append([key, fkey, family_size[fkey], family_score[fkey], rscore[key], dg[key], g4[key], structure_score[key], final_score[key], presentation[key], sequence[key]])
                nu += 1
        outf.close()
        df = pd.DataFrame(data, columns=["key", "fkey", "family_size", "family_score", "rscore", "dg", "g4", "structure_score", "final_score", "presentation", "sequence"])
        df.to_csv(self.output_dir + '/result_aptamer.csv')


if __name__ == '__main__':
    print("smart_mda_rnafold.py", str(datetime.datetime.now()))
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--kmer', type=int, default=6,
                        help='the predifined length of k-mers (default: 6)')
    parser.add_argument('-t', '--threads', type=int, default=multiprocessing.cpu_count(),
                        help='threads (default: device CPU number)')
    parser.add_argument('-c', '--temp', type=float, default=37,
                        help='rescale energy parameters to a temperature of temp C. (default  37)')
    parser.add_argument('-o', '--output', type=str, default="outputs",
                        help='output directory (default result)')
    args = parser.parse_args()
    k = args.kmer
    threads = args.threads
    temp = float(args.temp)
    output_dir = str(os.getcwd()) + "/" + args.output
    mdarnafolder = MDARNAFolder(k=k, threads=threads, temp=temp, output_dir=output_dir)
    mdarnafolder.run()