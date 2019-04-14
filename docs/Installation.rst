======================================
Installation
======================================

Required tools
-----------------
. topic:: Required tools
    * GATK >= 4.0 ( java 8 required )
    * Bedtools
    * Python 2.7 
    
        -scipy
    
        -numpy
    
        -multiprocessing
    
    * R
    
        -ztable
    
        -rmarkdown
    
        -corrplot


Installation with 'git clone'
------------------------------

::

    $ git clone https://github.com/heinc1010/BAMixChecker



Installation with the compressed file
---------------------------------------

Download the compressed file from https://github.com/heinc1010/BAMixChecker.

Then, decompress the .zip file.

::

    $ unzip BAMixChecker-master.zip


Set configeration
-------------------

And then set the tools PATH on the configuration file.


::

    $ cd BAMixChecker
    $ vim BAMixChecker.config
    
    GATK=/path/of/gatk
    BEDTOOLS=/path/of/bedtools

* Default (If you added the path in $PATH, you don’t need to modify the configuration file )
::

    GATK=gatk 
    BEDTOOLS=bedtools


