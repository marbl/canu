
.. _pipeline:

Canu Pipeline
=============

The canu pipeline is big and complicated and we haven't written it up yet.  Sorry.


The basic flow is:

- Gatekeeper
- Meryl
- Overlapper

  + Read Correction
  + Read Trimming

    * Trimming
    * Splitting

  + Unitig Construction

    * Overlap Error Adjustment
    * Overlap Filtering (NOT IMPLEMENTED)
    * Unitig
    * Consensus
    * Labeling
    * Output

Details:

Meryl - counts kmers, generates a histogram txt file, and a histogram plot png

