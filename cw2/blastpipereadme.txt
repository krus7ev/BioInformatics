blastpipe.py - Performs the blast pipeline processing

It requires the name of a FASTA <query file> to perform a search against. It also requires an ASTRAL (or other FASTA database) file to run, but it will also use the makeblastdb command to index the database if required.

Minimal Use:

py blastpipe.py -q <query file>

Will try to run against a database called "astral.1.75" in the current directory and then query it with <query file>.

Advanced use:

The pipeline can also use MUSCLE to create multiple sequence alignments and PhyML to generate phylogenetic trees from the results. When using this option you will need to use the "-m" and "-p" parameters (specified below). The program also requires the "fastatophyml.py" to be stored in the same directory as the program so that it can convert the FASTA files to the PhyML format.

Other options:
  -d <database name>: specifies the name and location of a database to use. The default is "astral.1.75". The .fa extension is ignored so you don't need to specify it (and it doesn't matter if you do). It looks for BLAST database files and if they don't exist it will try the makeblastdb command pointing at a .fa file with the same name.

  -m <muscle command>: Specify the path or command to execute muscle. On windows this will need to be the path to the muscle exe file. This parameter is required for the parts of the pipeline that perform sequence alignment.

  -p <phyml command>: Specify the path or command to execute PhyML. On Windows this will need to be the path to the exe file. This parameter is required for parts of the pipeline that create phylogenetic trees.

  -v: Verbose output. Will output the full names of the significant hits as results of the query.

  -i: Identifier output. Will output the identifiers of the significant hits as results of the query.

  -? <question number>: Performs the function indicated by the given question of the assignment. For example to run the program for question 17 send '-? 17'. The -i and -v arguments are overriden when you use the -? argument.
