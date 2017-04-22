
..
 introduction
  what this is
 
 quick start
  one PacBio SMRT cell from ecoli
 
 pipeline (not reference, designed to be read through)
  overview
   introduce modular pipeline
   introduce gkpStore, overlaps, ovlStore, tigStore
  read correction
  read trimming
  unitig construction
 
  local vs grid mode
 
 
 canu option reference
  each option, in detail, grouped by function
 
 canu executable reference
  each binary, in detail, alphabetical
 
 option index (alphabetical)

 history

Canu
====

.. toctree::
   :hidden:

   quick-start
   faq
   tutorial
   pipeline
   parameter-reference
   command-reference
   history


`Canu <http://github.com/marbl/canu>`_ is a fork of the Celera Assembler designed for high-noise single-molecule sequencing (such as
the PacBio RSII or Oxford Nanopore MinION). 

Publication
===========
Koren S, Walenz BP, Berlin K, Miller JR, Phillippy AM. `Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation <http://doi.org/10.1101/gr.215087.116>`_. Genome Research. (2017).

Install
=========
The easiest way to get started is to download a `release <https://github.com/marbl/canu/releases>`_. If you encounter
any issues, please report them using the `github issues <http://github.com/marbl/canu/issues>`_ page.

Alternatively, you can also build the latest unreleased from github:

::

  git clone https://github.com/marbl/canu.git
  cd canu/src
  make -j <number of threads>

Learn
=========

*  :ref:`Quick Start               <quickstart>` - no experience or data required, download and assemble *Escherichia coli* today!
*  :ref:`FAQ                       <faq>` - Frequently asked questions
*  :ref:`Canu tutorial             <tutorial>`   - a gentle introduction to the complexities of canu.
*  :ref:`Canu pipeline             <pipeline>`   - what, exactly, is canu doing, anyway?

*  :ref:`Canu Parameter Reference  <parameter-reference>` - all the paramters, grouped by function.
*  :ref:`Canu Command Reference    <command-reference>` - all the commands that canu runs for you.
*  :ref:`Canu History              <history>` - the history of the Canu project.
