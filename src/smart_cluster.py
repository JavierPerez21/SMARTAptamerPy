import argparse
import os
import multiprocessing
import subprocess
import datetime


class Cluster():
    def __init__(self, threads=multiprocessing.cpu_count(), output_dir=str(os.getcwd()) + "/outputs",
                 primer_file=str(os.getcwd()) + "/data/primers.txt", inflation=1.5, evalue=5e-2,
                 ):
        self.threads = threads
        self.output_dir = output_dir
        self.primer_file = primer_file
        self.inflation = inflation
        self.evalue = evalue

    def run(self):
        self.create_fvar_n_back()
        self.do_blast()
        self.gen_graph()
        self.perform_mcl()
        self.cluster_format_change()
        self.final_score_function()
        subprocess.call([f"rm -rf {self.output_dir + '/family_fasta'}"], shell=True)
        subprocess.call([f"rm -rf {self.output_dir + '/fastadir'}"], shell=True)


    def create_fvar_n_back(self):
        """
        Read primer file and create fvar and back files
        :return:
        """
        fvar = ''
        back = ''
        if os.path.exists(self.primer_file):
            with open(self.primer_file) as inf:
                nu = 0
                for line in inf:
                    if nu == 1:
                        back = line
                    else:
                        fvar = line
                        nu += 1
            inf.close()
        else:
            print("WARNING: This isn't any file including primers")
        if not os.path.exists(self.output_dir + "/family_fasta"):
            os.mkdir(self.output_dir + "/family_fasta")
        self.fvar = fvar.replace("\n", "")
        self.back = back.replace("\n", "")

    def do_blast(self):
        """
        Run blastn on fvar and back files
        :return:
        """
        print("***STEP: running blast***\n")
        subprocess.call([f"makeblastdb -in {self.output_dir}/all_used_uniq.fasta -dbtype nucl -out {self.output_dir}/goodsequence.fasta"], shell=True)
        subprocess.call([f"blastn -task blastn-short -evalue {self.evalue} -db {self.output_dir}/goodsequence.fasta -query {self.output_dir}/all_used_uniq.fasta -outfmt 7 -out {self.output_dir}/goodsequence_blastn.fasta -dust no -num_threads {self.threads}"], shell=True)
        subprocess.call([f"grep -v -P \"^#\" {self.output_dir}/goodsequence_blastn.fasta > {self.output_dir}/goodsequence_v1_blastn.out"], shell=True)

    def gen_graph(self):
        """
        Generate graph
        :return:
        """
        similiar1 = {}
        sim = {}
        with open(self.output_dir + "/goodsequence_v1_blastn.out") as inf:
            for line in inf:
                info = line.split("\t")
                if info[0] in similiar1 and info[1] in similiar1[info[0]]:
                    similiar1[info[0]][info[1]] += float(info[11])
                    similiar1[info[1]][info[0]] += float(info[11])
                    sim[info[0] + "-" + info[1]] += float(info[11])
                else:
                    similiar1[info[0]] = {info[1]: float(info[11])}
                    similiar1[info[1]] = {info[0]: float(info[11])}
                    sim[info[0] + "-" + info[1]] = float(info[11])
        inf.close()
        similiar2 = list(sim.values())
        maxim = max(similiar2)
        minim = min(similiar2)
        mean = sum(similiar2) / len(similiar2)
        self.mcl_map = {}
        self.r_mcl_map = {}
        n_similar = {}
        nu = 0
        with open(self.output_dir + "/all_used_uniq.fasta") as inf:
            for line in inf:
                if line.startswith(">"):
                    self.mcl_map[nu] = line.strip().split(">")[1]
                    self.r_mcl_map[line.strip().split(">")[1]] = nu
                    nu += 1
        inf.close()
        # Graph
        nu -= 1
        with open(self.output_dir + "/graph.txt", "w") as out:
            out.write("(mclheader\nmcltype matrix\ndimensions ")
            out.write(str(nu))
            out.write("x")
            out.write(str(nu))
            out.write("\n")
            out.write(")\n")
            out.write("(mclmatrix\nbegin\n\n")
            nu1 = 0
            while nu1 <= nu:
                out.write(str(nu1))
                keys = similiar1[self.mcl_map[nu1]].keys()
                number = len(keys)
                if number > 0:
                    for key in keys:
                        similiar1[self.mcl_map[nu1]][key] = similiar1[self.mcl_map[nu1]][key] / maxim
                        out.write("\t" + str(self.r_mcl_map[key]) + ":" + str(similiar1[self.mcl_map[nu1]][key]))
                out.write("\t$\n")
                nu1 += 1
            out.write(")\n")

    def perform_mcl(self):
        """
        Perform MCL
        :return:
        """
        print("***STEP: running mcl***\n")
        subprocess.call([f"mcl {self.output_dir}/graph.txt -I {self.inflation} -o {self.output_dir}/mcl_cluster -te {self.threads} -V all"], shell=True)

    def cluster_format_change(self):
        """
        Change cluster format
        :return:
        """
        print("***STEP: changing cluster format***\n")
        id = ""
        self.seq = {}
        with open(self.output_dir + "/all_used_uniq.fasta", "r") as inf:
            for line in inf:
                line = line.strip()
                if line.startswith(">"):
                    id = line[1:]
                else:
                    self.seq[id] = line
        inf.close()
        self.nfreq = {}
        with open(self.output_dir + "/all_used_info.txt", "r") as inf:
            for line in inf:
                line = line.strip()
                info = line.split("\t")
                self.nfreq[info[0]] = info[2]
        inf.close()
        self.freq = {}
        self.members = {}
        self.nmembers = {}
        with open(self.output_dir + "/mcl_cluster", "r") as inf:
            lines = "".join(inf.readlines())
            lines = lines.split("begin")[1]
            lines = ",".join(lines.split())
            lines = lines.split(",$,")
            lines = lines[:-1]  #c#
            for block in lines:
                if len(block) == 0:
                    continue
                info = block.split(",")
                if not info[0].isdigit():
                    continue
                fid = info[0]
                fid = "fid_" + fid
                #print('fid:', fid)  #c#
                #print(info)  #c#
                membersstring = None
                with open(self.output_dir + "/family_fasta/" + fid + ".fasta", "w") as outf:
                    for member in info:
                        member = self.mcl_map[int(member)]
                        outf.write(">" + member + "\n" + self.seq[member] + "\n")
                        if membersstring is not None:
                            membersstring += ":" + member
                        else:
                            membersstring = member
                outf.close()
                info = membersstring.split(":")
                un = len(info)
                self.nmembers[fid] = un
                self.members[fid] = membersstring
                for member in info:
                    if fid in self.freq.keys():
                        self.freq[fid] += int(self.nfreq[member])
                    else:
                        self.freq[fid] = int(self.nfreq[member])
        inf.close()

        with open(self.output_dir + "/aptamer_clusters", "w") as outf:
            for key in sorted(self.freq, key=self.freq.get, reverse=True):
                outf.write(key + "\t" + str(self.freq[key]) + "\t" + str(self.nmembers[key]) + "\t" + self.members[key] + "\n")
        outf.close()

    def max_value1(self, hash, keys):
        """
        Find the max value in the hash
        :param hash:
        :param keys:
        :return:
        """
        max = -10000000000000
        mkey = ""
        for key in keys:
            if hash[key] > max:
                max = hash[key]
                mkey = key
            elif hash[key] == max:
                mkey += ":" + key
        return mkey

    def max_value2(self, hash, keys):
        """
        Find the max value in the hash
        :param hash:
        :param keys:
        :return:
        """
        max = -10000000000000
        mkey = ""
        for key in keys:
            if hash[key] > max:
                max = hash[key]
                mkey = key
        return mkey

    def final_score_function(self):
        """

        :return:
        """
        f2seqs = {}
        with open(self.output_dir + "/aptamer_clusters", "r") as inf:
            for line in inf:
                line = line.strip()
                info = line.split("\t")
                f2seqs[info[0]] = info[3]
        inf.close()

        self.score = {}
        with open(self.output_dir + "/all_used_info.txt", "r") as inf:
            for line in inf:
                line = line.strip()
                info = line.split("\t")
                self.score[info[0]] = float(info[1])
                self.nfreq[info[0]] = int(info[2])
                self.seq[info[0]] = info[3]
        inf.close()
        k1 = len(f2seqs)
        k2 = list(f2seqs.keys())
        pl = int(k1 / self.threads)
        ll = k1 - pl * self.threads
        nu = 0
        nus = {}
        while nu < self.threads:
            mus = []
            nn = 0
            if nu < ll:
                while nn <= pl:
                    ln = nu * (pl + 1) + nn
                    mus.append(ln)
                    nn += 1
            elif nu == ll:
                while nn <= (pl - 1):
                    ln = nu * (pl + 1) + nn
                    mus.append(ln)
                    nn += 1
            else:
                while nn <= (pl - 1):
                    ln = ll * (pl + 1) + (nu - ll) * pl + nn
                    mus.append(ln)
                    nn += 1
            nus[nu] = mus
            nu += 1

        fastadir = self.output_dir + "/fastadir"
        if not os.path.exists(fastadir):
            os.mkdir(fastadir)
        pm = multiprocessing.Pool(self.threads)        # Not sure if actually using
        for nu in nus:
            mus = nus[nu]
            with open(fastadir + "/represent_" + str(nu) + ".fasta", "w") as outf:
                for ln in mus:
                    fid = k2[ln]
                    sids = f2seqs[fid]
                    sids = sids.split(":")
                    if len(sids) == 1:
                        id = fid + "_" + sids[0]
                        outf.write(id + "\t" + self.fvar + self.seq[sids[0]] + self.back + "\n")
                    else:
                        keys = self.max_value1(self.nfreq, sids)
                        keys = keys.split(":")
                        key = self.max_value2(self.score, keys)
                        id = fid + "_" + key
                        outf.write(id + "\t" + self.fvar + self.seq[key] + self.back + "\n")
                outf.close()
        pm.close()
        pm.join()

        subprocess.call([f"cat {fastadir}/represent_*.fasta >{self.output_dir}/represent.fasta"], shell=True)


if __name__ == '__main__':
    print("smart_cluster.py", str(datetime.datetime.now()))
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threads', type=int, default=multiprocessing.cpu_count(),
                        help='threads (default: device CPU number)')
    parser.add_argument('-o', '--output', type=str, default="outputs",
                        help='output directory (default result)')
    parser.add_argument('-i', '--inflation', type=float, default=1.5,
                        help='inflation value for mcl algorithm (default 1.5)')
    parser.add_argument('-e', '--evalue', type=float, default=0.05,
                        help='cutoff e-value after blast results (default 0.05)')
    parser.add_argument('-p', '--primer', type=str, default="data/primers_1.txt",
                        help='the primer file')
    args = parser.parse_args()
    threads = args.threads
    output_dir = str(os.getcwd()) + "/" + args.output
    inflation = float(args.inflation)
    evalue = float(args.evalue)
    primer_file = str(os.getcwd()) + "/" + args.primer

    cluster = Cluster(threads=threads, output_dir=output_dir, inflation=inflation, evalue=evalue, primer_file=primer_file)
    cluster.run()