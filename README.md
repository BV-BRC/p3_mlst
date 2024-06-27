# BV-BRC MLST Prediction

## Overview

This module contains the code for predicting MLST tags for input genomes. It is
invoked as part of the BV-BRC [genome annotation service](https://www.bv-brc.org/docs/quick_references/services/genome_annotation_service.html).

Multilocus sequence typing (MLST) is an unambiguous procedure for characterising isolates of bacterial species using the sequences of internal fragments of (usually) seven house-keeping genes. Approximately 450-500 bp internal fragments of each gene are used, as these can be accurately sequenced on both strands using an automated DNA sequencer. For each house-keeping gene, the different sequences present within a bacterial species are assigned as distinct alleles and, for each isolate, the alleles at each of the seven loci define the allelic profile or sequence type (ST). (Quoting from [PubMLST](https://pubmlst.org/multilocus-sequence-typing))

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

The script [`assign_st_to_genome`](scripts/assign_genomes_to_dbs.pl) is invoked during the execution of the BV-BRC genome annotation pipeline. The output of the script is parsed and used to assign a MLST sequence type for the genome being analyzed.

The libraries used are downloaded based on the index data from http://pubmlst.org/data/dbases.xml.