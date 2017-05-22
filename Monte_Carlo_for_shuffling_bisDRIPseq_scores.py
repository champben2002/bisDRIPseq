#Purpose: to generate a Monte Carlo simulation of our bisDRIP-seq experiments where the bisDRIP-seq scores of each read are decided by stochastic processes. The simulation works by randomly shuffling the location of bisDRIP-seq scores across the reads.

#Output of runswitchlabelsonchrms(): A multi-column tab-delimited file with bisDRIP-seq reads from a single chromosome associated with randomly assigned bisDRIP-seq scores. All reads in this file are associated with scores above zero. Each row of the file will contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row should contain the start position of each read. The third column in each row should contain the end position of each read. The final column of each row should contain the randomly-assigned bisDRIP-seq score. 

#run runswitchlabelsonchrms(readbisDRIPseqscorefile,allreadfile,outputfolder, name) using the following inputs:

#readbisDRIPseqscorefile:  A multi-column tab-delimited file with bisDRIP-seq reads from a single chromosome associated with bisDRIP-seq scores. Reads associated with bisDRIP-seq scores that equal zero may be excluded from this file. Each row of the file should contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row should contain the start position of each read. The third column in each row should contain the end position of each read. The final column of each row should contain the bisDRIP-seq score of the read.

#allreadfile:  A multi-column tab-delimited file with all bisDRIP-seq reads from a single chromosome. Each row of the file should contain information regarding a single read. Each row should have multiple columns.

#outputfolder: a directory where the output will be written into

#name - name of the file that the contents will be written into. 

#Outline of function:
#runswitchlabelsonchrms(readbisDRIPseqscorefile,allreadfile,outputfolder, name)
	#1 addlabelstoreads(rfile, cfile,workingof + "temp",workingof + "finaltemp.txt")
    	#Shuffles and adds shuffled bisDRIP-seq scores to reads
    		#1A filelength()
       		#determines how many reads aligned to genome.
	#2 reducesize(workingof + "finaltemp.txt", ofile)
    		#2A filelength()

import os
import random

def filelength(inputfile):
#1A Determines how many reads aligned to genome.
    a = open(inputfile, "r")
    b = a.readline()
    ic = 0
    c = 0
    t = ""
    while b != "" and ic < 100:
        e = b.replace("\n", "")
        f = e.split("\t")
        if len(f) < 2:
            ic += 1
        else:
            t = b
            c += 1
        b = a.readline()
    a.close()
    print(t)
    return(c)

def addlabelstoreads(scorefile, allreadfile,tempstem,outputfile):
#1 Shuffles and adds shuffled bisDRIP-seq scores to reads
    allreadlen = int(filelength(allreadfile))
    tempfilea = tempstem + "tempa.txt"
    tempfileb = tempstem + "tempb.txt"
    tempfilec = tempstem + "tempc.txt"
    tf = open(tempfilea, "w")
    a = open(scorefile, "r")
    b = a.readline()
    rset = []
    ic = 0
    first = 0
    while b != "" and ic < 100:
        e = b.replace("\n", "")
        f = e.split("\t")
        if len(f) < 2:
            ic += 1
        else:
            if first == 0:
                chk = 0
                while chk < 1:
                    rndnum = random.randint(0,allreadlen)
                    if rset.count(rndnum) == 0:
                        chk = 2
                rset.append(rndnum)
                tf.write(str(rndnum) + "\t" + f[-1])
                first = 1
            else:
                chk = 0
                while chk < 1:
                    rndnum = random.randint(0,allreadlen)
                    if rset.count(rndnum) == 0:
                        chk = 2
                rset.append(rndnum)
                tf.write("\n" + str(rndnum) + "\t" + f[-1])
        b = a.readline()
    a.close()
    tf.close()
    os.system("sort -k1,1n " + tempfilea + "> " + tempfileb)
    os.system("sort -k1,1 -k2,2n " + allreadfile + "> " + tempfilec)
    ff = open(outputfile, "w")
    a = open(tempfilec, "r")
    check = open(tempfileb, "r")
    cb = check.readline()
    ce = cb.replace("\n", "")
    cf = ce.split("\t")
    b = a.readline()
    ic = 0
    first = 0
    while b != "" and ic < 100:
        e = b.replace("\n", "")
        f = e.split("\t")
        if len(f) < 2:
            ic += 1
        else:
            if cf[0] == "":
                ff.write("\n" + "\t".join(f) + "\t0")
                print(cf)
            elif int(cf[0]) == first:
                print(first)
                if first == 0:
                    ff.write("\t".join(f) + "\t" + str(cf[1]))
                else:
                    ff.write("\n" + "\t".join(f) + "\t" + str(cf[1]))
                cb = check.readline()
                ce = cb.replace("\n", "")
                cf = ce.split("\t")
            else:
                if first == 0:
                    ff.write("\t".join(f) + "\t0")
                else:
                    ff.write("\n" + "\t".join(f) + "\t0")
            first += 1
        b = a.readline()
    a.close()
    ff.close()
    os.system("rm " + tempfilea)
    os.system("rm " + tempfileb)
    os.system("rm " + tempfilec)

def reducesize(inputfile,outputfile):
#2 Removes reads that have zero bisDRIP-seq scores
    a = open(inputfile, "r")
    ff = open(outputfile, "w")
    b = a.readline()
    ic = 0
    first = 0
    while b != "" and ic < 100:
        e = b.replace("\n", "")
        f = e.split("\t")
        if len(f) < 2:
            ic += 1
        else:
            if float(f[-1]) > 0:
                if first == 0:
                    ff.write(e)
                    first += 1
                else:
                    ff.write("\n" + e)
        b = a.readline()
    a.close()
    ff.close()

def runswitchlabelsonchrms(readbisDRIPseqscorefile,allreadfile,outputfolder, name):
    os.system("mkdir " + outputfolder)
    addlabelstoreads(readbisDRIPseqscorefile, allreadfile,outputfolder + "temp",outputfolder + "finaltemp.txt")
    reducesize(outputfolder + "finaltemp.txt", outputfolder + name)
    os.system("rm " + outputfolder + "finaltemp.txt")



