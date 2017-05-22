#Purpose: to generate a Monte Carlo simulation of our bisDRIP-seq experiments where the region containing each bisDRIP-seq read is determined by stochastic processes. The simulation works by randomly shuffling the location of bisDRIP-seq reads with non-zero scores.

#Output of runrandomizeronafileRRF(): A multi-column tab-delimited file with bisDRIP-seq reads from a single chromosome associated with bisDRIP-seq scores and a randomly assigned location. All reads in this file are associated with scores above zero. Each row of the file will contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row should contain the randomly-assigned start position of each read. The third column in each row should contain the assigned end position of each read. The final column of each row should contain the bisDRIP-seq score of the read. 

#run runrandomizeronafileRRF(readfile,checkfile,tempstem,outputfile) using the following inputs:

#readfile:  A multi-column tab-delimited file with bisDRIP-seq reads from a single chromosome associated with bisDRIP-seq scores. Reads associated with bisDRIP-seq scores that equal zero may be excluded from this file. Each row of the file should contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row should contain the start position of each read. The third column in each row should contain the end nucleotide of each read. The final column of each row should contain the bisDRIP-seq score of the read.

#checkfile:  A multi-column tab-delimited file with regions from a single chromosome that contain bisDRIP-seq reads. Each row of the file should contain information regarding a single region. The first column of each row must contain the chromosome of the region using the notation chrA where A is the number or letter of the chromosome. The second column of each region should contain the median positon of the region.

#tempstem: a directory and possible start of file name where files will be temporary written and read from

#outputfile - name (including directory) of the file that the contents will be written into. 

#Outline of function:
#runrandomizeronafileRRF(readfile,checkfile,tempstem,outputfile):
	#1 randomlyassignreads(readfile,checkfile,outputfile, "no")
   	#assign pseudoreads anywhere in a chromosome
	#2 absenteereads(outputfile,checkfile,tempstem,tempstem + "tempc.txt")
    	#Removes pseudoreads that do not fall in 1kb regions that reads align to
		#2A combinefeatureandregionfiles()
		#2B checkreadlistforregions()
	#3 nofilecontent(tempstem + "tempc.txt")

import os
import random
import statistics

def randomlyassignreads(realreadfile, evenreadfile,outputfile,alreadyopen):
#1 assign reads anywhere in a chromosome
    maxsize = 0
    a = open(evenreadfile, "r")
    b = a.readline()
    ic = 0
    while b != "" and ic < 100:
        e = b.replace("\n", "")
        f = e.split("\t")
        if len(f) < 2:
            ic += 1
        else:
            if float(f[1]) + 50 > maxsize:
                maxsize = float(f[1]) + 50
        b = a.readline()
    a.close()
    maxsize += 150
    if alreadyopen == "yes":
        ff = open(outputfile, "a")
    if alreadyopen == "no":
        ff = open(outputfile, "w")
    a = open(realreadfile, "r")
    b = a.readline()
    ic = 0
    while b != "" and ic < 100:
        e = b.replace("\n", "")
        f = e.split("\t")
        if len(f) < 2:
            ic += 1
        else:
            newstart = random.randint(0,maxsize)
            dif = int(f[2])-int(f[1])
            if alreadyopen == "no":
                ff.write(f[0] + "\t" + str(newstart) + "\t" + str(newstart + dif) + "\t" + "\t".join(f[3:]))
                alreadyopen = "yes"
            else:
                ff.write("\n" + f[0] + "\t" + str(newstart) + "\t" + str(newstart + dif) + "\t" + "\t".join(f[3:]))
        b = a.readline()
    a.close()
    ff.close()

def checkreadlistforregions(readlist, lowmax, regionmin):
#2B
    c = 0
    final = [[],[]]
    while c < len(readlist):
        if int(readlist[c][1]) > int(lowmax) and int(readlist[c][2]) < int(regionmin):
            final[1].append(readlist[c])
        else:
            final[0].append(readlist[c])
        c += 1
    return final

def combinefeatureandregionfiles(featurefile,regionfile, tempfile):
#2A
    a = open(featurefile, "r")
    ff = open(tempfile, "w")
    b = a.readline()
    b = a.readline()
    c = 0
    ic = 0
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            if c == 0:
                ff.write(f[0] + "\t" + "\t".join(f[1:]))
                c += 1
            else:
                ff.write("\n" + f[0] + "\t" + "\t".join(f[1:]))
        b = a.readline()
    a.close()
    a = open(regionfile, "r")
    b = a.readline()
    c = 0
    ic = 0
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            if float(f[2]) > 0 or float(f[3]) > 0:
                ff.write("\n" + "\t".join(f) + "\t" + "region")
        b = a.readline()
    a.close()
    ff.close()

def absenteereads(readfile,checkfile,tempstem,outputfile,distance):
#2
    tempa = tempstem + "tempa.txt"
    tempb = tempstem + "tempb.txt"
    combinefeatureandregionfiles(readfile,checkfile, tempa)
    os.system("sort -k1,1 -k2,2n " + tempa + "> " + tempb)
    a = open(tempb, "r")
    b = a.readline()
    first = 0
    firstt = 0
    ic = 0
    regionmin = 0
    lowmax = 0
    readlist = []
    final = []
    ff = open(outputfile, "w")
    fft = open(tempa, "w")
    while b != "" and ic < 100:
        e = b.replace("\n", "")
        f = e.split("\t")
        if len(f) < 2:
            ic += 1
        else:
            if f[-1] == "region":
                regionmin = int(f[1]) - 51
                final = checkreadlistforregions(readlist, lowmax, regionmin)
                c = 0
                while c < len(final[0]):
                    if firstt == 0:
                        fft.write("\t".join(final[0][c]))
                        firstt += 1
                    else:
                        fft.write("\n" + "\t".join(final[0][c]))
                    c += 1
                c = 0
                while c < len(final[1]):
                    if first == 0:
                        ff.write("\t".join(final[1][c]))
                        first += 1
                    else:
                        ff.write("\n" + "\t".join(final[1][c]))
                    c += 1
                lowmax = int(f[1]) + 50
                readlist = []
            else:
                readlist.append(f)
        b = a.readline()
        c = 0
    while c < len(readlist):
        if int(readlist[c][1]) > int(lowmax):
            ff.write("\n" + "\t".join(readlist[c]))
        else:
            fft.write("\n" + "\t".join(readlist[c]))
        c += 1
    a.close()
    ff.close()

def nofilecontent(inputfile):
    a = open(inputfile)
    empty = "false"
    b = a.readline()
    if b == "":
        b = a.readline()
        if b == "":
            empty = "true"
    return empty

def runrandomizeronafileRRF(readfile,checkfile,tempstem,outputfile):
    randomlyassignreads(readfile,checkfile,outputfile, "no")
    absenteereads(outputfile,checkfile,tempstem,tempstem + "tempc.txt")
    os.system("cp " + tempstem + "tempa.txt " + outputfile)
    while nofilecontent(tempstem + "tempc.txt") == "false":
        randomlyassignreads(tempstem + "tempc.txt",checkfile,tempstem + "tempa.txt", "yes")
        os.system("cp " + tempstem + "tempa.txt " + outputfile)
        absenteereads(outputfile,checkfile,tempstem,tempstem + "tempc.txt")
    os.system('rm ' + tempstem + "tempa.txt")
    os.system('rm ' + tempstem + "tempb.txt")
    os.system('rm ' + tempstem + "tempc.txt")
