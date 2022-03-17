import argparse
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import multiprocessing
import subprocess
import datetime
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF
import pandas as pd


class Scorer():
    def __init__(self, k=6, threads=multiprocessing.cpu_count(), quantile=0.995, control='theo06',
                 list_file=str(os.getcwd()) + '/data/theolist.txt',
                 input_dir=str(os.getcwd()) + '/data/theophylline_txt',
                 output_dir=str(os.getcwd()) + "/outputs",
                 directory=str(os.getcwd()) +"/src",
                 ecdf_per=0.9
                 ):
        self.k = k
        self.threads = threads
        self.quantile = quantile
        self.control = control
        self.list_file = list_file
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.directory = directory
        self.ecdf_per = ecdf_per

    def run(self):
        self.count_kmers_all_rounds()
        #self.run_R_kmer_count_pdf_script()
        self.get_kmer_count_pdf()
        self.create_count_enrich_kmers_dicts()
        self.generate_shift(self.count, self.rounds)
        self.kmer_score(self.shift, self.enrich_kmers)
        self.output_kmer_scores()

    def count_kmers_all_rounds(self):
        """
        Count kmer in each round listed in list_file
        :return:
        """
        print("Counting kmers...")
        self.rounds = []
        with open(self.list_file, 'r') as inf:
            for line in inf:
                self.rounds.append(line.strip())
        # Count kmer in each round
        p = multiprocessing.Pool(self.threads)
        for round in self.rounds:
            p.apply_async(self.count_kmers, args=(round,))
        p.close()
        p.join()
        print("Done counting kmers.")

    def count_kmers(self, round):
        """
        Count kmer in round
        :param round: round name as specified in list_file
        :return:
        """
        with open(self.input_dir + '/' + round + ".txt", 'r') as f:
            data = f.readlines()
        count = {}
        reads = {}
        for read in data:
            if read[0] == '>':
                continue
            else:
                seq = read.strip()
                lkmer = self.k
                nu = 0
                temp = {}
                while nu <= (len(seq) - lkmer):
                    substr = seq[nu:nu + lkmer]
                    if temp.get(substr):
                        pass
                    else:
                        if count.get(substr):
                            count[substr] += 1
                        else:
                            count[substr] = 1
                        temp[substr] = 1
                    nu += 1
                if reads.get(round):
                    reads[round] += 1
                else:
                    reads[round] = 1
        with open(self.output_dir + '/' + round + '_kmer.txt', 'w') as f:
            for kmer in count:
                f.write(kmer + '\t' + str(count[kmer]) + '\n')
        with open(self.output_dir + '/' + round + '_reads.txt', 'w') as f:
            f.write(str(reads[round]))

    def get_kmer_count_pdf(self):
        fix, axs = plt.subplots(int(len(self.rounds) / 2) + 1, 2, figsize=(10, 10))
        cut_offs = []
        for k, round in enumerate(self.rounds):
            read_file = self.output_dir + "/" + round + "_kmer.txt"
            data = pd.read_table(read_file, header=None)
            data[1] = data[1] / sum(data[1])
            axs[int(k / 2), k % 2].hist(data[1], bins=300)
            axs[int(k / 2), k % 2].set_title("Historgram of kmers")
            axs[int(k / 2), k % 2].set_xlabel(round)
            axs[int(k / 2), k % 2].set_ylabel("Frequency")
            with open(self.output_dir +"/" + round + "_kmer_percent.txt", "w") as outf:
                outf.write("\"V1\" \"V2\"\n")
                for i in range(0, len(data)):
                    row = f"\"{str(i + 1)}\" \"{data[0][i]}\" {data[1][i]}\n"
                    outf.write(row)
            outf.close()
            q95 = np.quantile(data[1], self.quantile)
            q95_info = round + "    " + str(q95)
            cut_offs.append(q95_info)
        # Adjust spacing between axs
        plt.tight_layout(pad=3.0, w_pad=1.0, h_pad=1.0)
        plt.savefig(self.output_dir + "/kmer.pdf")
        with open(self.output_dir +"/" + "cut_off.txt", "w") as outf:
            for i in range(0, len(cut_offs)):
                outf.write(cut_offs[i] + "\n")
        outf.close()

    def create_count_enrich_kmers_dicts(self):
        """
        Create dictionary of kmer counts for each round
        :return:
        """
        print("Creating count and enrich_kmers dictionaries...")
        with open(self.output_dir + '/cut_off.txt', 'r') as f:
            for line in f.readlines():
                if self.control in line:
                    #cutoff = float(line.split()[2][:-1])
                    cutoff = float(line.split()[1])
        self.count, self.enrich_kmers = {}, {}
        for round in self.rounds:
            with open(self.output_dir + '/' + round + '_kmer_percent.txt', 'r') as f:
                for line in f.readlines():
                    if "V1" in line and "V2" in line:
                        continue
                    info = line.strip().split()
                    if round == self.control:
                        if info[1] not in self.count.keys():
                            self.count[info[1]] = {round: float(info[2])}
                        else:
                            self.count[info[1]][round] = float(info[2])
                    else:
                        if float(info[2]) > cutoff:
                            self.enrich_kmers[info[1]] = 1
                        if info[1] not in self.count.keys():
                            self.count[info[1]] = {round: float(info[2])}
                        else:
                            self.count[info[1]][round] = float(info[2])
        print("Done creating count and enrich_kmers dictionaries.")

    def generate_shift(self, count, rounds):
        """
        Generate shift for each round
        :param count: self.count
        :param rounds: self.rounds
        :return:
        """
        rounds = [round.strip() for round in rounds]
        print("Generating shift...")
        length = len(rounds)
        length -= 1
        shift = {}
        for kmer in count:
            max = 0
            index = 0
            while index < length:
                ratio = 0
                if rounds[index] in count[kmer].keys():
                    if rounds[index+1] in count[kmer].keys():
                        ratio = count[kmer][rounds[index + 1]] / count[kmer][rounds[index]]
                    else:
                        ratio = 0
                else:
                    if rounds[index + 1] in count[kmer].keys():
                        ratio = 5
                    else:
                        ratio = 0
                if ratio > max:
                    max = ratio
                index += 1
            shift[kmer] = max
        self.shift = shift
        print("Done generating shift.")

    def kmer_score(self, shift, enrich_kmers):
        """
        Calculate score for each kmer
        :param shift: self.shift
        :param enrich_kmers: self.enrich_kmers
        :return:
        """
        print("Calculating kmer score...")
        shifts = list(shift.values())
        shifts.sort()
        nkmer = 4 ** self.k
        nkmer1 = len(enrich_kmers)
        nkmer = int(nkmer * (1 - nkmer1 / nkmer))
        print(nkmer, len(shifts), nkmer1)
        cutoff2 = shifts[nkmer]
        keys = list(shift.keys())
        for kmer in keys:
            if shift[kmer] > cutoff2:
                shift[kmer] = shift[kmer] - 1
            else:
                del shift[kmer]
        # Define score for each kmer
        score_kmers = {}
        for kmer in shift:
            score_kmers[kmer] = 1
        for kmer in self.enrich_kmers:
            score_kmers[kmer] = 1

        for kmer in score_kmers:
            if kmer in self.enrich_kmers.keys():
                score_kmers[kmer] += self.enrich_kmers[kmer]
            if kmer in self.shift.keys():
                score_kmers[kmer] += self.shift[kmer]

        self.score_kmers = score_kmers
        print("Done calculating kmer score.")

    def get_kmer_score_pdf(self):
        file = os.path.join(self.output_dir, "score_kmers.txt")
        data = pd.read_table(file, header=None)
        kmer_scores1 = data[1]
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fn1 = ECDF(kmer_scores1)
        x = fn1.x
        y = fn1.y
        self.filter_score_cutoff = x[np.where(y > self.ecdf_per)[0][0]]
        ax1.plot(fn1.x, fn1.y)
        ax1.scatter(x[np.where(y > self.ecdf_per)[0][0]], y[np.where(y > self.ecdf_per)[0][0]], c='r')
        ax1.annotate(f"({round(x[np.where(y > self.ecdf_per)[0][0]], 2)},{round(y[np.where(y > self.ecdf_per)[0][0]], 2)})",
                     (x[np.where(y > self.ecdf_per)[0][0]], y[np.where(y > self.ecdf_per)[0][0]]))
        ax1.set_title("Analyzed K-mer scores")
        ax1.set_xlabel("K-mer score")
        ax1.set_ylabel("Cumulative probability")
        ax1.grid(True)
        kmer_scores2 = np.append(kmer_scores1, np.zeros(4 ** self.k - len(kmer_scores1)))
        fn2 = ECDF(kmer_scores2)
        ax2.plot(fn2.x, fn2.y)
        ax2.set_title("All K-mer scores")
        ax2.set_xlabel("K-mer score")
        ax2.grid(True)
        plt.savefig(os.path.join(self.output_dir, "ecdf.pdf"))

    def output_kmer_scores(self):
        """
        Generate outputs of kmer scores
        :return:
        """
        print("Generating output of kmer scores...")
        # Generate score_kmers.txt
        with open(self.output_dir + "/score_kmers.txt", 'w') as f:
            for kmer in sorted(self.score_kmers, key=self.score_kmers.get, reverse=True):
                f.write(kmer[1:-1] + "\t" + str(self.score_kmers[kmer]) + "\n")
        self.get_kmer_score_pdf()
        print("Done generating output of kmer scores.")

if __name__ == '__main__':
    print("find_score.py", str(datetime.datetime.now()))
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--kmer', type=int, default=6,
                        help='the predifined length of k-mers (default: 6)')
    parser.add_argument('-t', '--threads', type=int, default=multiprocessing.cpu_count(),
                        help='threads (default: device CPU number)')
    parser.add_argument('-q', '--quantile', type=float, default=0.995,
                        help='the quantile that used to define the enriched k-mers (default: 0.995)')
    parser.add_argument('-c', '--control', type=str, default='theo06',
                        help='the control round (default R0)')
    parser.add_argument('-f', '--list', type=str, default='data/theolist_cut_16.txt',
                        help='library list (default list.txt)')
    parser.add_argument('-i', '--input', type=str, default='data/theophylline_txt_cut_16',
                        help='input directory (default input)')
    parser.add_argument('-o', '--output', type=str, default="outputs",
                        help='output directory (default result)')
    parser.add_argument('-d', '--directory', type=str, default="src",
                        help='source directory that contain the SMART-Apta')
    args = parser.parse_args()
    k = args.kmer
    threads = args.threads
    quantile = args.quantile
    control = args.control
    list_file = str(os.getcwd()) + "/" + args.list
    input_dir = str(os.getcwd()) + "/" + args.input
    output_dir = str(os.getcwd()) + "/" + args.output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    directory = str(os.getcwd()) + "/" + args.directory
    scorer = Scorer(k=k, threads=threads, quantile=quantile, control=control, list_file=list_file, input_dir=input_dir,
                    output_dir=output_dir, directory=directory)
    scorer.run()

