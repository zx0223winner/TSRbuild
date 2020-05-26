#Contributed by the AtPEAT study authors
#
matchAdapters

C++ program that will strip out the adapters.
Takes fasta (or fastq) format and outputs data in the same format, minus the
adapter. Only outputs read pairs that contain a 5' adapter on one read
and a 3' adapter on the other read.
