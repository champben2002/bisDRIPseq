import os

#Run allproximalscorestosinglentfeatures(featurefile,inputfolder,outputfolder,name, direction,disfromfeature) to get a metaplot of read bisDRIP-seq scores relative to a group of features in featurefile. The first column of the feature file must be the chromosome expressed as "chrX" where X is 1,2,3,...,X. The second column must be the location of the feature. Inputfolder must contain bisDRIPseq score files.The metaplot will calculate the bisDRIPseq score at each position relative to features. The distance is the limit to how far away from the feature that the metaplot will be calculated. Each feature cannot be closer than 2X distance from another feature. Output will be put into a tab-delimmited file in 'outputfolder + name + final.txt'. Chromosome will be in column 1, position relative to feature in 2 and bisDRIP-seq score in column 3 to final column. BisDRIP-seq scores are put into seperate columns for each inputfile.
#1 proximalscorestart()
#1A combinefeatureandreadfiles()
#1B saveproximalreads()
#1Bi quicklistthin()
#1C scoreregions()
#2 proximalscorecont()
#2A combinefeatureandreadfiles()
#2B saveproximalreads()
#2Bi quicklistthin()
#2C addscoreregions()

def combinefeatureandreadfiles(featurefile,readfile, tempfile):
#1A and 2A
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
#1Bi and 2Bi
    c = 0
    newlist = []
    while c < len(list):
        if int(minnt)-int(list[c][2]) < disfromfeature:
            newlist += [list[c]]
        c += 1
    return newlist

def saveproximalreads(tempfilea,tempfileb,disfromfeature):
#1B and 2B
    a = open(tempfilea, "r")
    ff = open(tempfileb, "w")
    b = a.readline()
    check = []
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
                    if int(f[1])-int(feature) < disfromfeature:
                        ff.write("\n" + "\t".join(f) + '\t' + str(int(f[1])-int(feature)) + '\t' + str(int(f[2])-int(feature)))
                    else:
                        ison = "no"
                        check = [f]
            else:
                c = 0
                if f[1] == "":
                    print(f)
                feature = int(f[1])
                st = quicklistthin(check,int(f[1]),disfromfeature)
                while c < len(st):
                    if firstrow == 0:
                        ff.write("\t".join(st[c]) + "\t" + str(int(st[c][1]) - feature)+ "\t" + str(int(st[c][2]) - feature))
                        firstrow += 1
                    else:
                        ff.write("\n" + "\t".join(st[c]) + "\t" + str(int(st[c][1]) - feature)+ "\t" + str(int(st[c][2]) - feature))
                    c += 1
                if firstrow == 0:
                    ff.write("\t".join(f))
                    firstrow += 1
                else:
                    ff.write("\n" + "\t".join(f))
                ison = "yes"
                chr = f[0]
        b = a.readline()
    a.close()
    ff.close()

def scoreregions(tempfilec,finalfile, filen, direction,disfromfeature):
#1C
    scorelist = []
    c = -disfromfeature
    while c < disfromfeature+1:
        scorelist.append([c,0])
        c += 1
    a = open(tempfilec, "r")
    ff = open(finalfile, "w")
    b = a.readline()
    ic = 0
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            if f[-3] == "read":
                c = 0
                while c < len(scorelist):
                    if direction == "positive":
                        if int(f[-2]) <= scorelist[c][0] <= int(f[-1]):
                            scorelist[c][1] += float(f[-4])
                    elif direction == "negative":
                        if int(f[-2]) <= -1*scorelist[c][0] <= int(f[-1]):
                            scorelist[c][1] += float(f[-4])
                    else:
                        print("ERROR with direction")
                    c += 1
        b = a.readline()
    ff.write("Relative position\t" + filen)
    c = 0
    while c < len(scorelist):
        ff.write("\n" + str(scorelist[c][0]) + "\t" + str(scorelist[c][1]))
        c += 1
    a.close()
    ff.close()

def proximalscorestart(featurefile,inputfile, folder, name, direction,disfromfeature):
#1
    tempa = folder + name + "ta.txt"
    tempb = folder + name + "tb.txt"
    tempc = folder + name + "tc.txt"
    ff = folder + name + "final.txt"
    combinefeatureandreadfiles(featurefile,inputfile, tempa)
    os.system("sort -k1,1 -k2,2n " + tempa + "> " + tempb)
    saveproximalreads(tempb,tempc, disfromfeature)
    os.system("rm " + tempa)
    os.system("rm " + tempb)
    scoreregions(tempc, ff, inputfile,direction, disfromfeature)
    os.system("rm " + tempc)

def addscoreregions(tempfilec,tempfiled, finalfile,filen, direction,disfromfeature):
#2C
    os.system("cp " + finalfile + " " + tempfiled)
    scorelist = []
    c = -disfromfeature
    while c < disfromfeature+1:
        scorelist.append([c,0])
        c += 1
    a = open(tempfilec, "r")
    ff = open(finalfile, "w")
    b = a.readline()
    ic = 0
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            if f[-3] == "read":
                c = 0
                while c < len(scorelist):
                    if direction == "positive":
                        if int(f[-2]) <= scorelist[c][0] <= int(f[-1]):
                            scorelist[c][1] += float(f[-4])
                    elif direction == "negative":
                        if int(f[-2]) <= -1*scorelist[c][0] <= int(f[-1]):
                            scorelist[c][1] += float(f[-4])
                    else:
                        print("ERROR with direction")
                    c += 1
        b = a.readline()
    a = open(tempfiled, "r")
    b = a.readline()
    e = b.replace("\n", "")
    ff.write(e + "\t" + filen)
    b = a.readline()
    ic = 0
    c = -disfromfeature
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            print(b + "\t" + str(scorelist[int(f[0])][0]) + "\t" + str(scorelist[int(f[0])][1]))
            ff.write("\n" + e + "\t" + str(scorelist[c + disfromfeature][1]))
            c += 1
        b = a.readline()
    a.close()
    ff.close()

def proximalscorecont(featurefile,inputfile, folder, name, direction,disfromfeature):
#2
    tempa = folder + name + "ta.txt"
    tempb = folder + name + "tb.txt"
    tempc = folder + name + "tc.txt"
    tempd = folder + name + "td.txt"
    ff = folder + name + "final.txt"
    combinefeatureandreadfiles(featurefile,inputfile, tempa)
    os.system("sort -k1,1 -k2,2n " + tempa + "> " + tempb)
    saveproximalreads(tempb,tempc,disfromfeature)
    os.system("rm " + tempa)
    os.system("rm " + tempb)
    addscoreregions(tempc,tempd, ff, direction + inputfile,direction,disfromfeature)
    os.system("rm " + tempc)
    os.system("rm " + tempd)

def allproximalscorestosinglentfeatures(featurefile,inputfolder,outputfolder,name, direction,disfromfeature):
    os.system("mkdir " + outputfolder)
    start = os.listdir(inputfolder)
    start.sort()
    c = 0
    proximalscorestart(featurefile,inputfolder + start[c], outputfolder, name, direction,disfromfeature)
    c = 1
    while c < len(start):
        proximalscorecont(featurefile,inputfolder + start[c], outputfolder, name, direction,disfromfeature)
        c += 1
        print(c)
    print("finished")
