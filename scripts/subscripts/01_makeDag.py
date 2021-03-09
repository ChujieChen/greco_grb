#!/usr/bin/python
"""
Make the dag for 01.sub
"""
import sys
import os

import glob

####### Check before running ########
MODE = "production"           # production or testing
script_path = '/scratch/cjchen/submit/01.sub'
## [::-1] reverse order
tws_in_second    = [10,       25,    50,   100,   250,   500][::-1]   
njobs            = [0, 0, 0, 0, 0, 10000][::-1]
batchsizes       = [1000,   1000,  1000,  1000,  1000,  1000][::-1]
grb_name         = "GRB180423A"    # an example grb used to create controled seeds
ncpu             = 2       # multi-core may result in out-of-memory
#######################################

dagname = '01.dag'
dag_contents = ""

for (tw_in_second, nbatch, batchNtrials) in zip(tws_in_second, njobs, batchsizes):
    totalNtrials = batchNtrials * nbatch
    mode = MODE
    for batchIndex in range(nbatch):
        jobname = "{}_batchSize{}_batchIndex{}_tw{}".format(grb_name, batchNtrials, batchIndex, tw_in_second)
        # In condor, every option needs quotes around it
        dag_contents += "JOB %s %s\n" % (jobname, script_path)
        dag_contents += "VARS %s " % (jobname)
        dag_contents += " grb_name=\"{}\"".format(grb_name)
        dag_contents += " batchNtrials=\"{}\"".format(batchNtrials)
        dag_contents += " batchIndex=\"{}\"".format(batchIndex)
        dag_contents += " tw_in_second=\"{}\"".format(tw_in_second)
        dag_contents += " totalNtrials=\"{}\"".format(totalNtrials)
        dag_contents += " ncpu=\"{}\"".format(ncpu)
        dag_contents += " mode=\"{}\"".format(mode)
        dag_contents += " outfilename=\"{}_batchSize{}_batchIndex{}_tw{}.npz\"".format(grb_name, 
                                                                    batchNtrials, 
                                                                    batchIndex, 
                                                                    tw_in_second)
        dag_contents += "\n"
       
               
if dag_contents != "":
    dagman = open('./'+dagname, "w")
    dagman.write(dag_contents)
    dagman.close()
    
####################################
#
# Don't forget to copy .sub and .dag 
# to /scratch/cjchen/submit/ before run:                     
#     condor_q_dag 01.dag    
#
####################################