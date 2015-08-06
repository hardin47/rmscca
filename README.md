rmscca
======

robust multiple sparce canonical correlation analysis

The functions needed in the analysis are the following:


All programs
------------

* scca.CVperm (does the cross validation and permutation -- all the work to get over curse of dimensionality)
* scca.function (does the thresholding for calculating lambda values.  Very important function for scca.CVperm)
* sample.sigma12.function (calculates the cross covariance matrix, K)



Only simulating (both null and correlated data)
-----------------------------------------------

* Cov.suped  (makes the correlation matrices for X and Y)
* sim.setup (simulates the data!)
* build.B  (generates the relationships between X and Y)


Interpreting results
--------------------

* parse.breast (count things from the breast cancer results)
* interpret.results.curveonly  (uses the Q-permutation curve to count positives / negatives)
* results (more parsing output)
* results.helper (more parsing output)
* determine.true.vals (for the blocking)



###  To run  ###

# Null: nullSimScript.R
	- generates null data with k=0 (i.e., B==0).  Only uses the first CC to determine if anything is called significant.

# Breast cancer data:  dataScript2.R (only chrom2)  dataScript.R (all chrom in parallel)
	- Uses the .csv files from the PMA package.  
	- the main function doing the work is scca.CVperm
	- the output contains the list of coefficients and correlations for both spearman and pearson

# Simulated data: fullSimScript.R
	- I need to work through all the parsing scripts again.  But the basic idea is to simulate data, run the RMSCCA code, then count things like "complete groups", false positives, etc.




