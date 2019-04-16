======================================
Issues
======================================

Running time of BAMixChecker
-------------------------------

In our paper, we reported the running time of BAMixChecker which showed 5.3 min with 4 processes and 9.9 min with 1 process for WES. 

The same non-linear speed tendancy also showed in records for the targeted sequencing dataset.

We believe the non-linear speed is due to CPU resource limitation in the tested environment. 

To apply on a controlled environment without other jobs, we tested using a desktop with Intel® Core™ i7-4790 CPU 3.60GHz with quad cores and 32 GB memory. 

The problem is GATK runs fast (approximately 20 sec per WES sample on the environment) but on the start of the program, CPU usage get high for a second and goes down after that. 

With the desktop, we observed CPU usage is peaked to 58% with 1 process already (not only by GATK, additional to basic system running) and it goes down to 20-30%. 

So compete for CPU between processes occurs even with 2 processes in the environment at peak. 

The total advantage could not applied. During the time using normal range of CPU, each process also were competed with 4 processes. 

Even though we couldn’t test in a best environment which we can utilize the advantage of BAMixChecker in, we reported the record because BAMixChecker still showed enough fast running time. 


+---------------------+-----+-----+-----+-----+
| # of Process	      |1	  |2	  |3	  |4    |
+---------------------+-----+-----+-----+-----+
| Running time(min) 	|9.9	|6.4	|5.7	|5.3  |
+---------------------------------------------+


Table R1. Average running time of BAMixChecker with various number of processes for 10 times with 30 WES samples in a desktop with Intel® Core™ i7-4790 CPU 3.60GHz with quad cores and 32 GB memory.

We recommand the user to adjust the number of process ccording to their environment.
