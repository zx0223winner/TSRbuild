This directory contains all instructions for reproducing the data work
described in

	Ganote, C. & Brendel, V.P (2019)
	Identification and characterization of plant promoters
	The Plant Cell, to be submitted

As preliminary steps, please follow the instructions in

	ZmGENOME/0README

After the preliminary steps are taken care of, all the data work can be
reproduced by executing

```bash
./xrun-with-simg
```

which will download the TSRchitect Singularity image and execute the workflow
scripts within the container. What exactly will be run is determined by flags
in the xrun-with-simg file. By default, pulltsrsimg=1. The other flags
correspond to the various data sets analyzed and are set to 0 (see the 0README
files in the subdirectories for documentation). Probably best would be to run
one experiment at a time, but the script will happily do everything if you ask
it to do so.

NOTE: You should definitely take a look at the scripts and configuration files
      to customize the workflow for your computing environment and needs. For
      example, ZmCAGE/ZmCAGE.configfile specifies

	THREADS=8

which you should change to the number of processors you wish to use.

NOTE: You will need about 60G of disk space to run the ZmCAGE workflow.

Enjoy

Contact:	Volker Brendel (vbrendel@indiana.edu)
