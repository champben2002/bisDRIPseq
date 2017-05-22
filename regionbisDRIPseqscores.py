#Purpose: to calculate bisDRIP-seq score for each region provided in a file ie the region from the transcription start site and +1000 bp for HIST1H1E, RPPH1, and MALAT1.

# Run the file proximalscorestart() to get the bisDRIP-seq score for each region as described below:

#Output of proximalscorestart(): The output is a multi-column tab-delimited file. Each row contains a different feature of interest. In the first column of each row is the chromosome of the feature. In the second column of each row is the position of the feature. In the third column of each row is the bisDRIP-seq score associated with that feature.


#Input of proximalscorestart(featurefile,inputfile, folder, name, direction,disfromfeature):

#featurefile - a tab-delimited file that contains the location of each feature. The first row must be a header. Each following row must contain a feature and each feature must be on a seperate row. The first column of each feature row must contain the chromosome of the feature using the notation chrA where A is the number or letter of the chromosome. The second column of each feature row must contain the genomic position of the feature.

#inputfile - a bisDRIP-seqs read score file. The contents of the file should contain bisDRIP-seq reads with their associated bisDRIP-seq scores. Each row of the files should contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row should contain the start position of each read. The third column in each row should contain the end position of each read. The final column of each row should contain the bisDRIP-seq score. 

#folder - directory where output should be written

#name - output will be written into the file: folder + name + "final.txt"

#Direction: if "positive" than the output provides the score after and including the feature otherwise you get back the score before and including the feature

#disfromfeature: The returned score will be of the sum of the bisDRIP-seq scores in the region between the feature and the disfromfeature


#Outline of functions:
#proximalscorestart(featurefile,inputfile, folder, name, direction,disfromfeature)
#1 combinefeatureandreadfiles(featurefile,inputfile, tempa)
    #combines the input read score file with the feature file and adds a column containing the word "read" to then end of the input read score file
#2 os.system("sort -k1,1 -k2,2n " + tempa + "> " + tempb)
    #Sorts the file by chromosome then start site
#3 saveproximalreads(tempb,ff,disfromfeature,direction)
    #Calculates the region score for each feature in the feature file
    #3A quicklistthin()
        #removes reads and read scores from consideration if they are not proximal to the next feature
    #3B sumgenescores
        #Sums scores of reads near the studied feature, when is then written to the outputfile

import os

def combinefeatureandreadfiles(featurefile,readfile, tempfile):
#combines the input read score file with the feature file and adds a column containing the word "read" to then end of the input read score file
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
    a = open(readfile, "r")
    b = a.readline()
    ic = 0
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            ff.write("\n" + "\t".join(f) + "\t" + "read")
        b = a.readline()
    a.close()
    ff.close()

def quicklistthin(list,minnt, disfromfeature):
#3A - removes reads and read scores from consideration if they are not proximal to the next feature
    c = 0
    newlist = []
    while c < len(list):
        if int(minnt)-int(list[c][2]) < disfromfeature:
            newlist += [list[c]]
        c += 1
    return newlist

def sumgenescores(genelist,feature,disfromfeature):
#3B - Sums scores of reads near the studied feature, when is then written to the outputfile
    c = 0
    score = [0,0]
    featurerangepos = range(feature,feature + disfromfeature+1)
    featurerangeneg = range(feature-disfromfeature,feature+1)
    while c < len(genelist):
        readrange = range(int(genelist[c][1]),int(genelist[c][2])+1)
        rset = set(readrange)
        print(genelist[c])
        score[0] += float(genelist[c][-4])*len(rset.intersection(featurerangeneg))/(1+int(genelist[c][2])-int(genelist[c][1]))
        score[1] += float(genelist[c][-4])*len(rset.intersection(featurerangepos))/(1+int(genelist[c][2])-int(genelist[c][1]))
        c += 1
    return score

def saveproximalreads(tempfilea,tempfileb,disfromfeature,direction):
#Calculates the region score for each feature in the feature file
    a = open(tempfilea, "r")
    ff = open(tempfileb, "w")
    b = a.readline()
    check = []
    randcount = 0
    chr = "chr"
    ison = "no"
    firstrow = 0
    ic = 0
    feature = -200
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            if f[-1] == "read":
                if ison == "no":
                    if chr != f[0]:
                        check = [f]
                        chr = f[0]
                    if chr == f[0]:
                        check = quicklistthin(check, int(f[1]),disfromfeature) + [f]
                if ison == "yes":
                    if int(f[1])-int(feature) < disfromfeature and chr == f[0]:
                        check += [f + [str(int(f[1])-int(feature)),str(int(f[2])-int(feature))]]
                    else:
                        ison = "no"
                        if firstrow == 0:
                            sumgenessc = sumgenescores(check,feature,disfromfeature)
                            if direction == "positive":
                                ff.write(str(chr) + "\t" + str(feature) + "\t" + str(sumgenessc[1]))
                            if direction == "negative":
                                ff.write(str(chr) + "\t" + str(feature) + "\t" + str(sumgenessc[0]))
                            firstrow += 1
                        else:
                            sumgenessc = sumgenescores(check,feature,disfromfeature)
                            if direction == "positive":
                                ff.write("\n" + str(chr) + "\t" + str(feature) + "\t" + str(sumgenessc[1]))
                            if direction == "negative":
                                ff.write("\n" + str(chr) + "\t" + str(feature) + "\t" + str(sumgenessc[0]))
                        check = [f]
            else:
                skip = 0
                if ison == "yes":
                    if f[0] == chr and str(f[1])==str(feature):
                        skip = 1
                    else:
                        sumgenessc = sumgenescores(check,feature,disfromfeature)
                        if direction == "positive":
                            if firstrow == 0:
                                ff.write(str(chr) + "\t" + str(feature) + "\t" + str(sumgenessc[1]))
                                firstrow += 1
                            else:
                                ff.write("\n" + str(chr) + "\t" + str(feature) + "\t" + str(sumgenessc[1]))
                        if direction == "negative":
                            if firstrow == 0:
                                ff.write(str(chr) + "\t" + str(feature) + "\t" + str(sumgenessc[0]))
                                firstrow += 1
                            else:
                                ff.write("\n" + str(chr) + "\t" + str(feature) + "\t" + str(sumgenessc[0]))
                        ncount = 0
                        while ncount < len(check):
                            check[ncount] = check[ncount][:-2]
                            ncount += 1
                if skip == 0:
                    c = 0
                    feature = int(f[1])
                    st = quicklistthin(check,int(f[1]),disfromfeature)
                    check = []
                    if f[0] == chr:
                        while c < len(st):
                            check.append(st[c] + [str(int(st[c][1]) - feature), str(int(st[c][2]) - feature)])
                            c += 1
                    else:
                        check = []
                        chr = f[0]
                ison = "yes"
        b = a.readline()
        randcount += 1
    sumgenessc = sumgenescores(check,feature,disfromfeature)
    if direction == "positive":
        ff.write("\n" + str(chr) + "\t" + str(feature) + "\t" + str(sumgenessc[1]))
    if direction == "negative":
        ff.write("\n" + str(chr) + "\t" + str(feature) + "\t" + str(sumgenessc[0]))
    a.close()
    ff.close()

def proximalscorestart(featurefile,inputfile, folder, name, direction,disfromfeature):
    print(featurefile)
    print(inputfile)
    print(name)
    tempa = folder + name + "ta.txt"
    tempb = folder + name + "tb.txt"
    ff = folder + name + "final.txt"
    combinefeatureandreadfiles(featurefile,inputfile, tempa)
    print("combine complete")
    os.system("sort -k1,1 -k2,2n " + tempa + "> " + tempb)
    print("sort complete")
    saveproximalreads(tempb,ff,disfromfeature,direction)
    print("reads saved")
    os.system("rm " + tempa)
    os.system("rm " + tempb)

proximalscorestart(featurefile,inputfile, folder, name, direction,disfromfeature)
