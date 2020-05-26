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
* [TSRbuild](https://github.com/BrendelGroup/TSRbuild) (this repository)
* [MoVRs](https://github.com/BrendelGroup/MoVRs)

Here we walk you through the data work described in our publications:

[Policastro2020](./Policastro2020)
Robert A. Policastro, R. Taylor Raborn, Volker P. Brendel, and Gabriel E. Zentner
(2019) _Simple and efficient mapping of transcription start sites with STRIPE-seq._
In review.

[Ganote2020](./Ganote2020)
Ganote, C. & Brendel, V.P
(2020) _Identification and characterization of plant promoters._
The Plant Cell, to be submitted

**TSRexploreR**
Robert A. Policastro _et al._
(2020) TSRexploreR: _A tool set for evaluation and characterization of transcription start regions in eukaryotic genomes._
In preparation.

Please see instructions on how to proceed in the named subdirectories.
