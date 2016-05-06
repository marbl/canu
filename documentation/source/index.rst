
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
the PacBio RSII or Oxford Nanopore MinION). You can `download <http://github.com/marbl/canu/releases>`_ a release. If you encounter
any issues, please report them using the `github issues <http://github.com/marbl/canu/issues>`_ page.

*  :ref:`Quick Start               <quickstart>` - no experience or data required, download and assemble *Escherichia coli* today!
*  :ref:`FAQ                       <faq>` - Frequently asked questions
*  :ref:`Canu tutorial             <tutorial>`   - a gentle introduction to the complexities of canu.
*  :ref:`Canu pipeline             <pipeline>`   - what, exactly, is canu doing, anyway?

*  :ref:`Canu Parameter Reference  <parameter-reference>` - all the paramters, grouped by function.
*  :ref:`Canu Command Reference    <command-reference>` - all the commands that canu runs for you.
*  :ref:`Canu History              <history>` - the history of the Canu project.
