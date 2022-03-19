# SMARTAptamerPy

Python implementation of [SMART Aptamer](https://github.com/JavierPerez21/SMARTAptamerPy.git) from [this paper](https://pubs.acs.org/doi/abs/10.1021/acs.analchem.9b05203).

This will allow others to use it to make modifications to it easier.



## Dependencies

The following programs are required to run this program:

* [RNAfold](https://algosb2019.sciencesconf.org/data/RNAtutorial.pdf)
* [MCL](http://micans.org/mcl/)
* [QGRS](https://github.com/freezer333/qgrs-cpp)

## Instructions

1. Install the dependencies listed above.
2. Install libraries in requirements.txt
3. Put fastq files for each round in data/targ_fastq_?
4. Write name of each round file in data/targ_fastq_? in seperate lines without .fastq
5. Write each primer in a separate line in data/targ_primers.txt
6. Run fastqcleaner.py
7. Run fastq2txt.py
8. Run main.py