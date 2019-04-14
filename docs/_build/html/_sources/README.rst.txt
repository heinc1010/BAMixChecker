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

It’s simple and fast but accurate to detect a pair of bam files from same individual with WGS/WES, RNA, Targeted sequencing dataset. 

And the user can catch the information of mismatched sample information as well as the matched sample information at a glance.

.. image:: Report_ex.gif
    :alt: BAMixChecker.html
    :width: 600

.. image:: Heatmap_ex.jpg
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




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
