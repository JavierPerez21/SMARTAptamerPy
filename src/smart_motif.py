import argparse
import os
import multiprocessing
import subprocess
import datetime


class MotifFinder():
    def __init__(self, k=6, threads=multiprocessing.cpu_count(), list_file=str(os.getcwd()) + '/data/theolist.txt',
                 t_content=0.6, score_cutoff=10, seq_num=15000,
                 input_dir=str(os.getcwd()) + '/data/theophylline_txt',
                 output_dir=str(os.getcwd()) + "/outputs",
                 directory=str(os.getcwd()) + "/src"
                 ):
        self.k = k
        self.threads = threads
        self.t_content = t_content
        self.score_cutoff = score_cutoff
        self.seq_num = seq_num
        self.list_file = list_file
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.directory = directory
        with open(self.list_file, 'r') as f:
            self.rounds = f.readlines()
        self.rounds = [round.strip() for round in self.rounds]

    def run(self):
        self.count_kmers_2()
        self.get_polyt()
        subprocess.call([f"rm -rf {self.output_dir+'/seeds'}"], shell=True)

    def score_func(self, seq, lkmer, nu=0, score=0):
        scores = []
        while nu < (len(seq) - lkmer):
            kmer = seq[nu:nu + lkmer]
            if kmer in self.score_kmers.keys() and self.score_kmers[kmer] > score:
                score = self.score_kmers[kmer]
            nu += 1
        return score

    def count_kmers(self):
        print("Counting kmers")
        self.score_kmers = {}
        with open(self.output_dir + '/score_kmers.txt', 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                self.score_kmers[line[0]] = float(line[1])

        with open(self.output_dir + '/all.txt', 'w') as outf:
            for round in self.rounds:
                with open(self.input_dir + '/' + round + '.txt', 'r') as inf:
                    for line in inf:
                        outf.write(line)


        seedsdir = self.output_dir + '/seeds'
        if not os.path.exists(seedsdir):
            os.makedirs(seedsdir)
        subprocess.call(['split', '-l', str(self.threads), '-d', '--additional-suffix=seed', self.output_dir + '/all.txt', seedsdir+"/"])

        tseeds = [seed for seed in os.listdir(seedsdir) if seed.endswith('seed')]
        for file in tseeds:
            newinfile = file.strip()
            newinfile = newinfile.split('_')
            newoutfile = seedsdir + "/" + newinfile[0] + '_score.txt'
            newinfile = seedsdir + "/" +newinfile[0] + '_uniq.txt'
            file = seedsdir + "/" + file
            subprocess.call([f'cat {file}|sort|uniq > {newinfile}'], shell=True)
            with open(newinfile, 'r') as inf:
                with open(newoutfile, 'w') as outf:
                    for line in inf:
                        seq = line.strip()
                        seq = seq.split('\t')
                        outf.write(seq[0] + '\t' + str(self.score_func(seq[0], self.k)) + '\n')

            subprocess.call([f"cat {seedsdir}/*seed_score.txt > {self.output_dir}/scores.txt"], shell=True)
        print("Finished counting kmers")

    def count_kmers_2(self):
        print("Counting kmers")
        self.score_kmers = {}
        with open(self.output_dir + '/score_kmers.txt', 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                self.score_kmers[line[0]] = float(line[1])

        with open(self.output_dir + '/all.txt', 'w') as outf:
            for round in self.rounds:
                with open(self.input_dir + '/' + round + '.txt', 'r') as inf:
                    for line in inf:
                        outf.write(line)

        with open(self.output_dir + '/all.txt', 'r') as inf:
            with open(self.output_dir + '/scores.txt', 'w') as outf:
                for line in inf:
                    seq = line.strip()
                    seq = seq.split('\t')
                    outf.write(seq[0] + '\t' + str(self.score_func(seq[0], self.k)) + '\n')
        outf.close()
        inf.close()
        print("Finished counting kmers")

    def get_polyt(self):
        """
        Gets the proportion of the genome that is poly-T
        :return:
        """
        print("Getting poly-T proportion")
        score = {}
        Tfilter = 0
        freq = {}
        with open(self.output_dir + '/scores.txt', 'r') as inf:
            with open(self.output_dir + '/polyT.txt', 'w') as outf:
                for line in inf:
                    line = line.strip()
                    info = line.split('\t')
                    seq = info[0]
                    countT = seq.count('T')
                    countT = countT / len(seq)
                    if countT >= self.t_content:
                        outf.write(line + '\n')
                        Tfilter += 1
                        continue
                    score[seq] = info[1]
        inf.close()
        outf.close()
        print(Tfilter)

        for round in self.rounds:
            with open(self.input_dir + '/' + round + '.txt', 'r') as inf:
                for line in inf:
                    line = line.strip()
                    if line in score.keys():
                        if line in freq.keys():
                            freq[line] += 1
                        else:
                            freq[line] = 1
        inf.close()

        nu = 0
        with open(self.output_dir + '/all_used_uniq.fasta', 'w') as outf:
            with open(self.output_dir + '/all_used_info.txt', 'w') as outf2:
                if self.score_cutoff is not None:
                    keys = sorted(score.keys(), key=lambda x: float(score[x]), reverse=True)
                    for seq in keys:
                        if float(score[seq]) > self.score_cutoff:
                            outf.write('>seq_' + str(nu) + '\n' + seq + '\n')
                            outf2.write('seq_' + str(nu) + '\t' + score[seq] + '\t' + str(freq[seq]) + '\t' + seq + '\n')
                            nu += 1
                            if nu > self.seq_num:
                                break
                elif self.seq_num is not None:
                    keys = sorted(score.keys(), key=lambda x: float(score[x]), reverse=True)
                    score_cutoff = float(score[keys[self.seq_num]])
                    for seq in score.keys():
                        if float(score[seq]) >= score_cutoff:
                            outf.write('>seq_' + str(nu) + '\n' + seq + '\n')
                            outf2.write('seq_' + str(nu) + '\t' + score[seq] + '\t' + str(freq[seq]) + '\t' + seq + '\n')
                            nu += 1
        outf.close()
        outf2.close()
        print("Finished getting poly-T proportion")


if __name__ == '__main__':
    print("smart_motif.py", str(datetime.datetime.now()))
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--kmer', type=int, default=6,
                        help='the predifined length of k-mers (default: 6)')
    parser.add_argument('-t', '--threads', type=int, default=multiprocessing.cpu_count(),
                        help='threads (default: device CPU number)')
    parser.add_argument('-f', '--list', type=str, default='data/theolist_cut_16.txt',
                        help='library list (default list.txt)')
    parser.add_argument('-p', '--t_content', type=float, default=0.6,
                        help='sequences with T bases > this cutoff are considered as t-rich sequences (0-1, default 0.6)')
    parser.add_argument('-s', '--score_cutoff', type=int, default=6,
                        help='the filter scores cutoff value used to filter sequences based on motif enrichment status (only need to set one between -s and -n, default 10)')
    parser.add_argument('-n', '--seq_num', type=int, default=15000,
                        help='the total unique sequences should be left after filtering (set one between -s and -n, default 15000)')
    parser.add_argument('-i', '--input', type=str, default='data/theophylline_txt_cut_16',
                        help='input directory (default input)')
    parser.add_argument('-o', '--output', type=str, default="outputs",
                        help='output directory (default result)')
    parser.add_argument('-d', '--directory', type=str, default="src",
                        help='source directory that contain the SMART-Apta')
    args = parser.parse_args()
    k = args.kmer
    threads = args.threads
    list_file = str(os.getcwd()) + "/" + args.list
    t_content = args.t_content
    score_cutoff = args.score_cutoff
    seq_num = args.seq_num
    input_dir = str(os.getcwd()) + "/" + args.input
    output_dir = str(os.getcwd()) + "/" + args.output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    directory = str(os.getcwd()) + "/" + args.directory
    motiffinder = MotifFinder(k=k, threads=threads, list_file=list_file, t_content=t_content, score_cutoff=score_cutoff,
                              seq_num=seq_num, input_dir=input_dir, output_dir=output_dir, directory=directory)
    motiffinder.run()