# Ready to work

Exploration of TSR presupposes that you have TSR to explore.
There are different experimental techniques for 5'-end transcript profiling.
Our computational analyses should be applicable to sequencing data from any
such techniques (both single- and paired-end read data).
The workflow involves the following basic steps:

* mapping of the read data to the relevant genome
* identifying TSS and TSR
* exploring the TSR
* (possibly) identification of transcription regulatory motifs

The following github repositories provide tools and documentation for all of
the steps (in order):

* [GoSTRIPES](https://github.com/BrendelGroup/GoSTRIPES)
* [TSRchitect](https://github.com/BrendelGroup/TSRchitect) (also available as a [Bioconductor](http://bioconductor.org/) package)
* [TSRexplore](https://github.com/BrendelGroup/TSRexplore) (this repository)
* [MoVRs](https://github.com/BrendelGroup/MoVRs)

Here we walk you through the data work described in our publications:

Robert A. Policastro, R. Taylor Raborn, Volker P. Brendel, and Gabriel E. Zentner
(2019) _Simple and efficient mapping of transcription start sites with STRIPE-seq._
To be submitted.

Robert A. Policastro _et al._
(2019) TSRexplore: _A tool set for evaluation and characterization of transcription start regions in eukaryotic genomes._
In preparation.

Give yourself a clean working directory on a nice Linux machine (10 cores,
100GB RAM, 300GB disk space recommended) and follow the instructions, starting
with

```bash
./xgetstarted1
```

All the files and the directory structure in the _work_ directory that contains
this _README_ file are necessary.
