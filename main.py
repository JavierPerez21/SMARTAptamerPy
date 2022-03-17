from src.smart_cluster import *
from src.find_score import *
from src.smart_motif import *
from src.smart_mda_rnafold import *


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--kmer', type=int, default=6,
                        help='the predifined length of k-mers (default: 6)')
    parser.add_argument('-t', '--threads', type=int, default=multiprocessing.cpu_count(),
                        help='threads (default: device CPU number)')
    parser.add_argument('-q', '--quantile', type=float, default=0.995,
                        help='the quantile that used to define the enriched k-mers (default: 0.995)')
    parser.add_argument('-c', '--control', type=str, default='theo06',
                        help='the control round (default R0)')
    parser.add_argument('-f', '--list', type=str, default='data/theolist_cut_18.txt',
                        help='library list (default list.txt)')
    parser.add_argument('-i', '--input', type=str, default='data/theophylline_txt',
                        help='input directory (default input)')
    parser.add_argument('-o', '--output', type=str, default="outputs",
                        help='output directory (default result)')
    parser.add_argument('-d', '--directory', type=str, default="src",
                        help='source directory that contain the SMART-Apta')
    parser.add_argument('-p', '--t_content', type=float, default=0.6,
                        help='sequences with T bases > this cutoff are considered as t-rich sequences (0-1, default 0.6)')
    parser.add_argument('-n', '--seq_num', type=int, default=50000,
                        help='the total unique sequences should be left after filtering (set one between -s and -n, default 15000)')
    parser.add_argument('-r', '--primer', type=str, default="data/primers_1.txt",
                        help='the primer file')
    parser.add_argument('-l', '--inflation', type=float, default=1.5,
                        help='inflation value for mcl algorithm (default 1.5)')
    parser.add_argument('-e', '--evalue', type=float, default=0.05,
                        help='cutoff e-value after blast results (default 0.05)')
    parser.add_argument('-m', '--temp', type=float, default=37,
                        help='rescale energy parameters to a temperature of temp C. (default  37)')
    parser.add_argument('-s', '--ecdf_per', type=float, default=0.9,
                        help='ECDF to cutoff to choose kmer_score_cutoff')

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
    t_content = args.t_content
    seq_num = args.seq_num
    primer_file = str(os.getcwd()) + "/" + args.primer
    inflation = float(args.inflation)
    evalue = float(args.evalue)
    temp = float(args.temp)
    ecdf_per = float(args.ecdf_per)

    with open(output_dir + "/config.txt", "w") as outf:
        outf.write("kmer: " + str(k) + "\n")
        outf.write("threads: " + str(threads) + "\n")
        outf.write("quantile: " + str(quantile) + "\n")
        outf.write("control: " + str(control) + "\n")
        outf.write("list: " + str(list_file) + "\n")
        outf.write("input: " + str(input_dir) + "\n")
        outf.write("output: " + str(output_dir) + "\n")
        outf.write("directory: " + str(directory) + "\n")
        outf.write("t_content: " + str(t_content) + "\n")
        outf.write("seq_num: " + str(seq_num) + "\n")
        outf.write("primer: " + str(primer_file) + "\n")
        outf.write("inflation: " + str(inflation) + "\n")
        outf.write("evalue: " + str(evalue) + "\n")
        outf.write("temp: " + str(temp) + "\n")
        outf.write("ecdf_per: " + str(ecdf_per) + "\n")
    outf.close()

    find_score_time = datetime.datetime.now()
    print("Running find_score.py...")
    scorer = Scorer(k=k, threads=threads, quantile=quantile, control=control, list_file=list_file, input_dir=input_dir,
                    output_dir=output_dir, directory=directory, ecdf_per=ecdf_per)
    scorer.run()

    smart_motif_time = datetime.datetime.now()
    print("find_score.py finished in " + str(smart_motif_time - find_score_time))
    print("Running smart_motif.py...")
    motiffinder = MotifFinder(k=k, threads=threads, list_file=list_file, t_content=t_content, score_cutoff=scorer.filter_score_cutoff,
                              seq_num=seq_num, input_dir=input_dir, output_dir=output_dir, directory=directory)
    motiffinder.run()

    smart_cluster_time = datetime.datetime.now()
    print("smart_motif.py finished in " + str(smart_cluster_time - smart_motif_time))
    print("Running smart_cluster.py...")
    cluster = Cluster(threads=threads, output_dir=output_dir, inflation=inflation, evalue=evalue,
                      primer_file=primer_file)
    cluster.run()

    smart_mda_rnafolder_time = datetime.datetime.now()
    print("smart_cluster.py finished in " + str(smart_mda_rnafolder_time - smart_cluster_time))
    print("Running smart_mda_rnafolder.py...")
    mdarnafolder = MDARNAFolder(k=k, threads=threads, temp=temp, output_dir=output_dir)
    mdarnafolder.run()

    final_time = datetime.datetime.now()
    print("smart_mda_rnafolder.py finished in " + str(final_time - smart_mda_rnafolder_time))
    print("All finished in " + str(final_time - find_score_time))

