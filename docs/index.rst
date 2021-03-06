.. BAMixChecker documentation master file, created by
   sphinx-quickstart on Sun Apr 14 15:45:58 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

================================================================================
BAMixChecker
================================================================================

: An automated tool for sample matching checkup in NGS cohorts
---------------------------------------------------------------

.. image:: Workflow.jpg
    :alt: Workflow of BAMixChecker
    :width: 700

BAMixChecker is a fast and easy-to-use sample matching checkup tool for NGS dataset.

It is simple and fast but accurately detects pairs of WGS, WES, RNA, targeted sequencing BAM/CRAM files originating from the same individual.

It informs the user about matched or mismatched sample at a glance.

**Report of BAMixChecker in HTML file**

.. image:: Report_ex.gif
    :alt: BAMixChecker.html
    :width: 600

**Heatmap result of BAMixChecker in pdf file**

.. image:: Heatmap_ex.gif
    :alt: BAMixChecker_heatmap.pdf
    :width: 500

.. toctree::
   :maxdepth: 3
   :caption: BAMixChecker:
   
   Installation
       Required tools
       Installation with 'git clone'
       Installation with the compressed file
       Set configeration

   Input
       BAM file
           RNA-seq BAM file
       Reference sequence
       SNP list for non-human organism
    
   Usage
       Required arguments
       Optional arguments
       Usage for each data type
       Tutorial
   
   Output
       Report of BAMixChecker
       Heatmap
       TEXT files
   
   Issues
       Running time of BAMixChecker
   
   Citation
   
   Contact


