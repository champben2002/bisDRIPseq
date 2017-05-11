import os

#calculate normalized bisDRIP-seq scores and cytosine normalized bisDRIP-seq scores by running runallalignment(). The folder of fastq read files is called ngsfold. The Normalized bisDRIP-seq scores of each read is placed in outputfolder + name + "CTOBstdnormalized.txt" and outputfolder + name + "CTOTstdnormalized.txt". Cytosine Normalized bisDRIP-seq scores of each read is placed in outputfolder + "cnorm" + name + "CTOBstdnormalized.txt" and outputfolder + "cnorm" + name + "CTOTstdnormalized.txt". Raw bisDRIP-seq scores of each read is placed in outputfolder + name + "CTOBwithscores.txt" and outputfolder + name + "CTOBwithscores.txt". Conversions at each file are placed in outputfolder + name + "CTOBconversions.txt" and outputfolder + name + "CTOTconversions.txt".

#runbisDRIPseqpipeline(ngsfold, outputfolder, name):
#1-5 aligns reads to genome
#   (1) runflexbar(ngsfold)
#       (1A)intflexbar()
#   (2) runbismark(ngsfold)
#        (2A) intbismark()
#   (3) prepareunmappedfull(ngsfold)
#   (4) secondrunbismark(ngsfold)
#       (4A) secondintbismark()
#   (5) finalrunbismark(ngsfold)
#       (5A) prepareunmappedfinal()
#       (5B) runfinalflexbar()
#           (5Bi) finalflexbar()
#       (5C) finalintbismark()
#6-9 determines which cytosines were converted
#   (6) runmethylationextractor(ngsfold)
#       (6A) preparemethylationfolder()
#   (7) rundepulication(ngsfold)
#   (8) actualextract(ngsfold)
#   (9) executedatadump(ngsfold)
#       (9A) wholefolderdatadump()
#           (9Ai) chrfolderdatadump()
#               (9Ai1) filedump()
#10-11 calculates raw bisDRIP-seq scores for each read
#   (10) combinereadsinfolderx(ngsfold)
#       (10A) combinereadsinfilex()
#   inputfile = ngsfold + "bismarkfiles/bismarkmethyloutput/chroutput/CTOBreads.txt"
#   (11) addbisDRIPseqscorestofile(inputfile, outputfolder, name + "CTOB")
#       (11A) conversionsbyindcsinallrloops()
#       (11B) calculaterandomconversion()
#           (11Bi) calculateavgconversion()
#       (11C) inversefile()
#       (11D) readscore()
#           (11Di) getinputdata()
#           (11Dii) expectedtablebinomial()
#               (11Dii1) binomial()
#                   (11Dii1a) choose()
#           (11Diii) rowreadscore()
#           (11Div) combinelists()
#       (11E) addscorestoinitialfile()
#12 and 13 create normalized bisDRIP-seq scores for reads using either standard normalization (12) or cytosine normalization (13)
#   (12) normreadscores(outputfolder + name + "CTOBwithscores.txt",outputfolder + name + "CTOBstdnormalized.txt")
#   (13) runallnormbyC(inputfile, outputfolder, "cnorm" + name + "CTOB")
#       (13A) normreadscorebyCs()
#           (13Ai) normscoresbyCs()
#   (11repeat) addbisDRIPseqscorestofile(inputfile, outputfolder, name + "CTOT")
#   (12repeat) normreadscores(outputfolder + name + "CTOTwithscores.txt",outputfolder + name + "CTOTstdnormalized.txt")
#   (13repeat) runallnormbyC(inputfile, outputfolder, "cnorm" + name + "CTOT")
#14 generates files that contain each cytosine position, the number of observed conversions of that cytosine, the fraction of cytosines converted, and the total number of reads containing that conversion.
#   (14) makeconversionsfile(ngsfold + "bismarkfiles/bismarkmethyloutput/chroutput/CTOB/", outputfolder + name + "CTOBconversions.txt")
#       (14A) addtobed
#           (14Ai) wbedline
#   (14repeat) makeconversionsfile(ngsfold + "bismarkfiles/bismarkmethyloutput/chroutput/CTOT/", outputfolder + name + "CTOTconversions.txt")

def intflexbar(ngsfilea,ngsfileb, ngsfolder, anum):
#1A
    os.system('/opt/flexbar_v2.5_linux64/flexbar -r ' + ngsfolder + ngsfilea + ' -p ' + ngsfolder + ngsfileb + ' -f i1.8 --pre-trim-phred 30 --pre-trim-left 10 -n 4 -t ' + ngsfolder + 'flexfiles/' + ngsfileb[:-10].lower()+anum)


def runflexbar(ngsfolder):
#1
    folderlist = os.listdir(ngsfolder)
    ngsfolderb = ngsfolder + 'flexfiles/'
    os.system('mkdir ' + ngsfolderb)
    counter = 0
    clist = []
    while counter < len(folderlist):
        if folderlist[counter][-1] == 'q':
            anum = int(folderlist[counter][-9:-6])
            if clist.count(anum) == 0:
                intflexbar(folderlist[counter][:-11]+'1'+folderlist[counter][-10:],folderlist[counter][:-11]+'2'+folderlist[counter][-10:],ngsfolder,folderlist[counter][-9:-6])
                counter += 1
                clist.append(anum)
            else:
                counter += 1
        else:
            counter += 1

def intbismark(ngsfilea,ngsfileb, ngsfold, anum):
#2A
    ngsfolder = ngsfold + 'flexfiles/'
    ngsfolderb = ngsfold + 'bismarkfiles/'
    os.system('/home/jason/bismark_v0.14.3/bismark --bowtie2 --pbat --samtools_path /opt/samtools-1.2.0/bin/ -X 1000 --unmapped -N 0 --multicore 2 --path_to_bowtie /home/jason/bowtie2-2.2.5/ -o ' + ngsfolderb + ' /home/jason/Documents/chroms/ -1 ' + ngsfolder + ngsfilea + ' -2 ' + ngsfolder + ngsfileb)

def runbismark(ngsfold):
#2
    ngsfolder = ngsfold + 'flexfiles/'
    folderlist = os.listdir(ngsfolder)
    os.system('mkdir ' + ngsfold + 'bismarkfiles/')
    counter = 0
    clist = []
    while counter < len(folderlist):
        if folderlist[counter][-1] == 'q':
            anum = int(folderlist[counter][-11:-8])
            if clist.count(anum) == 0:
                intbismark(folderlist[counter][:-7]+'1'+folderlist[counter][-6:],folderlist[counter][:-7]+'2'+folderlist[counter][-6:],ngsfold,folderlist[counter][-11:-8])
                counter += 1
                clist.append(anum)
            else:
                counter += 1
        else:
            counter += 1

def prepareunmappedfull(ngsfold):
#3
    os.system('mkdir ' + ngsfold + 'bismarkfiles/unmappedfull/')
    os.system('mkdir ' + ngsfold + 'bismarkfiles/bismarkfilesunmappedfull/')
    os.system('mv ' + ngsfold + 'bismarkfiles/*_unmapped_reads_* ' + ngsfold + 'bismarkfiles/unmappedfull/')
    os.system('gunzip ' + ngsfold + 'bismarkfiles/unmappedfull/*')

def secondintbismark(filea, ngsfolderb):
#4A
    os.system('/home/jason/bismark_v0.14.3/bismark --bowtie2 --pbat --samtools_path /opt/samtools-1.2.0/bin/ --unmapped -N 0 --multicore 2 --path_to_bowtie /home/jason/bowtie2-2.2.5/ -o ' + ngsfolderb + ' /home/jason/Documents/chroms/ ' + filea)

def secondrunbismark(ngsfold):
#4
    ngsfolder = ngsfold + 'bismarkfiles/unmappedfull/'
    ngsfolderb = ngsfold + 'bismarkfiles/bismarkfilesunmappedfull/'
    folderlist = os.listdir(ngsfolder)
    counter = 0
    while counter < len(folderlist):
        secondintbismark(ngsfolder + folderlist[counter], ngsfolderb)
        counter += 1

def prepareunmappedfinal(ngsfold):
#5A
    os.system('mkdir ' + ngsfold + 'bismarkfiles/unmapped55/')
    os.system('mv ' + ngsfold + 'bismarkfiles/bismarkfilesunmappedfull/*_unmapped_reads.fq* ' + ngsfold + 'bismarkfiles/unmapped55/')
    os.system('gunzip ' + ngsfold + 'bismarkfiles/unmapped55/*')

def finalflexbar(ngsfile, ngsfold, ngsfolderb):
#5Bi
    os.system('/opt/flexbar_v2.5_linux64/flexbar -r ' + ngsfold + 'bismarkfiles/unmapped55/' + ngsfile + ' -f i1.8 --pre-trim-phred 30 -k 55 -n 10 -t ' + ngsfolderb + ngsfile[:-3] + 'unmapped_55')

def runfinalflexbar(ngsfold):
#5B
    folderlist = os.listdir(ngsfold + 'bismarkfiles/unmapped55/')
    ngsfolderb = ngsfold + 'flexfilesunmappedfinal/'
    os.system('mkdir ' + ngsfolderb)
    counter = 0
    while counter < len(folderlist):
        if folderlist[counter][-1]== 'q':
            finalflexbar(folderlist[counter], ngsfold, ngsfolderb)
            counter += 1
        else:
            counter += 1

def finalintbismark(fileset, ngsfolder):
#5C
    os.system('/home/jason/bismark_v0.14.3/bismark --bowtie2 --pbat --samtools_path /opt/samtools-1.2.0/bin/ --unmapped -N 0 --multicore 2 --path_to_bowtie /home/jason/bowtie2-2.2.5/ -o ' + ngsfolder + ' /home/jason/Documents/chroms/ ' + fileset)

def finalrunbismark(ngsfold):
#5
    prepareunmappedfinal(ngsfold)
    runfinalflexbar(ngsfold)
    ngsfolder = ngsfold + 'flexfilesunmappedfinal/'
    folderlist = os.listdir(ngsfolder)
    os.system('mkdir ' + ngsfold + 'bismarkfiles/unmapped55bismarkfiles/')
    counter = 0
    while counter < len(folderlist):
        finalintbismark(ngsfolder + folderlist[counter], ngsfold + 'bismarkfiles/unmapped55bismarkfiles/')
        counter += 1

def preparemethylationfolder(ngsfold):
#6A
    methylpeinput = ngsfold + 'bismarkfiles/bismarkbamoutput/peinput/'
    methylsingleinput = ngsfold + 'bismarkfiles/bismarkbamoutput/singleinput/'
    methyloutput = ngsfold + 'bismarkfiles/bismarkmethyloutput/'
    os.system('mkdir ' + ngsfold + 'bismarkfiles/bismarkbamoutput/')
    os.system('mkdir ' + methylpeinput)
    os.system('mkdir ' + methylsingleinput)
    os.system('mkdir ' + methyloutput)
    os.system('mv ' + ngsfold + 'bismarkfiles/*.bam ' + methylpeinput)
    os.system('mv ' + ngsfold + 'bismarkfiles/bismarkfilesunmappedfull/*.bam ' + methylsingleinput)
    os.system('mv ' + ngsfold + 'bismarkfiles/unmapped55bismarkfiles/*.bam ' + methylsingleinput)

def runmethylationextractor(ngsfold):
#6
    preparemethylationfolder(ngsfold)
    methylpeinput = ngsfold + 'bismarkfiles/bismarkbamoutput/peinput/'
    methylsingleinput = ngsfold + 'bismarkfiles/bismarkbamoutput/singleinput/'
    methyloutput = ngsfold + 'bismarkfiles/bismarkmethyloutput/'
    methylpeinputlist = os.listdir(methylpeinput)
    methylsingleinputlist = os.listdir(methylsingleinput)
    counterpe = 1
    filesetpe = methylpeinput + methylpeinputlist[0]
    while counterpe < len(methylpeinputlist):
        filesetpe += ' ' + methylpeinput + methylpeinputlist[counterpe]
        counterpe += 1
    countersin = 1
    filesetsin = methylsingleinput + methylsingleinputlist[0]
    while countersin < len(methylsingleinputlist):
        filesetsin += ' ' + methylsingleinput + methylsingleinputlist[countersin]
        countersin += 1
    os.system('/home/jason/bismark_v0.14.3/bismark_methylation_extractor -o ' + methyloutput + ' -p --no_header --no_overlap --multicore 3 --samtools_path /opt/samtools-1.2.0/bin/ --merge_non_CpG ' + filesetpe)
    os.system('/home/jason/bismark_v0.14.3/bismark_methylation_extractor -o ' + methyloutput + ' -s --no_header --multicore 3 --samtools_path /opt/samtools-1.2.0/bin/ --ignore_3prime 10 --merge_non_CpG ' + filesetsin)

def rundepulication(ngsfold):
#7
    ngsfolder = ngsfold + "bismarkfiles/bismarkbamoutput/peinput/"
    counter = 0
    a = os.listdir(ngsfolder)
    while counter < len(a):
        job = ('/home/jason/bismark_v0.14.3/deduplicate_bismark --bam --samtools_path /opt/samtools-1.2.0/bin/ -p ' + ngsfolder + a[counter])
        os.system(job)
        counter += 1
    ngsfolder = ngsfold + "bismarkfiles/bismarkbamoutput/singleinput/"
    counter = 0
    a = os.listdir(ngsfolder)
    while counter < len(a):
        job = ('/home/jason/bismark_v0.14.3/deduplicate_bismark --bam --samtools_path /opt/samtools-1.2.0/bin/ -s ' + ngsfolder + a[counter])
        os.system(job)
        counter += 1

def actualextract(ngsfold):
#8
    methylpeinput = ngsfold + 'bismarkfiles/bismarkbamoutput/peinput/'
    methylsingleinput = ngsfold + 'bismarkfiles/bismarkbamoutput/singleinput/'
    methyloutput = ngsfold + 'bismarkfiles/bismarkmethyloutput/'
    withduplicate = ngsfold + 'bismarkfiles/withduplicate/'
    os.system("mkdir " + withduplicate)
    os.system("mv " + methylpeinput + "* " + withduplicate)
    os.system("mv " + withduplicate + "*duplicated* " + methylpeinput)
    os.system("mv " + methylsingleinput + "* " + withduplicate)
    os.system("mv " + withduplicate + "*duplicated* " + methylsingleinput)
    methylpeinputlist = os.listdir(methylpeinput)
    methylsingleinputlist = os.listdir(methylsingleinput)
    counterpe = 1
    filesetpe = methylpeinput + methylpeinputlist[0]
    while counterpe < len(methylpeinputlist):
        filesetpe += ' ' + methylpeinput + methylpeinputlist[counterpe]
        counterpe += 1
    countersin = 1
    filesetsin = methylsingleinput + methylsingleinputlist[0]
    while countersin < len(methylsingleinputlist):
        filesetsin += ' ' + methylsingleinput + methylsingleinputlist[countersin]
        countersin += 1
    methyloutput = ngsfold + 'bismarkfiles/bismarkmethyloutput/'
    os.system('/home/jason/bismark_v0.14.3/bismark_methylation_extractor -o ' + methyloutput + ' -p --no_header --no_overlap --multicore 3 --samtools_path /opt/samtools-1.2.0/bin/ --merge_non_CpG ' + filesetpe)
    os.system('/home/jason/bismark_v0.14.3/bismark_methylation_extractor -o ' + methyloutput + ' -s --no_header --multicore 3 --samtools_path /opt/samtools-1.2.0/bin/ --ignore_3prime 10 --merge_non_CpG ' + filesetsin)

def filedump(chrfolder,bfile):
#9Ai1
    check = os.listdir(chrfolder)
    wecounter = 0
    tsfile = open(bfile, 'r')
    ts = "start"
    while ts != "" or wecounter < 1000:
        ts = tsfile.readline()
        tsfile2 = ts.replace('\n','')
        tsfile3 = tsfile2.replace('"','')
        fl = tsfile3.split('\t')
        if len(fl) < 4:
            wecounter += 1
        else:
            if fl[2] + ".txt" in check:
                temper = open(chrfolder + fl[2] + ".txt", "a")
                temper.write("\n" + fl[3] + "\t" + fl[1] + "\t" + fl[4] + "\t" + fl[0])
                temper.close()

            else:
                temper = open(chrfolder + fl[2] + ".txt", "w")
                temper.write(fl[3] + "\t" + fl[1] + "\t" + fl[4] + "\t" + fl[0])
                temper.close()
                check.append(fl[2] + ".txt")
    tsfile.close()

def chrfolderdatadump(chrfolder, listoffiles):
#9Ai
    counter = 0
    while counter < len(listoffiles):
        filedump(chrfolder,listoffiles[counter])
        counter += 1

def wholefolderdatadump(chrfolder,filefolder):
#9A
    nlist = []
    Blist = []
    start = os.listdir(filefolder)
    counter = 0
    while counter < len(start):
        if start[counter][-20:] == 'bt2.deduplicated.txt':
            if start[counter].count("CTOT") == 1:
                nlist.append(filefolder + start[counter])
                counter += 1
            elif start[counter].count("CTOB") == 1:
                Blist.append(filefolder + start[counter])
                counter += 1
            else:
                counter += 1
        elif start[counter][-20:] == '_pe.deduplicated.txt':
            if start[counter].count("CTOT") == 1:
                nlist.append(filefolder + start[counter])
                counter += 1
            elif start[counter].count("CTOB") == 1:
                Blist.append(filefolder + start[counter])
                counter += 1
            else:
                counter += 1
        else:
            counter += 1
    chrfolderdatadump(chrfolder + 'CTOT/', nlist)
    chrfolderdatadump(chrfolder + 'CTOB/', Blist)

def executedatadump(ngsfold):
#9 Dump info into files for each chrom
    print(ngsfold)
    inputfolder = ngsfold + "bismarkfiles/bismarkmethyloutput/"
    os.system('mkdir ' + inputfolder + 'chroutput/')
    os.system('mkdir ' + inputfolder + 'chroutput/CTOT/')
    os.system('mkdir ' + inputfolder + 'chroutput/CTOB/')
    wholefolderdatadump(inputfolder + "chroutput/", inputfolder)

def combinereadsinfilex(filex, finalfile,readcol,chromosome):
#10A
    os.system("sort -k"+str(readcol+1)+","+str(readcol+1)+" -k1,1n -f " + filex + " -T /home/jason/ -o " + filex)
    a = open(filex,"r")
    b = open(finalfile,"a")
    d = a.readline()
    c = 0
    ci = 0
    check = ["",0,0,0,0,0,0,""]
    while d != "" and ci < 100:
        e = d.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ci += 1
        else:
            if check[7] == "":
                print(check)
                print(f)
                print(c)
                if f[1] == "+":
                    check = [chromosome,f[0],f[0],"-1","-1","0","1",f[readcol]]
                else:
                    check = [chromosome,f[0],f[0],f[0],f[0],"1","1",f[readcol]]
            elif f[readcol] == check[7]:
                if f[1] == "+":
                    check[2] = f[0]
                    check[6] = str(int(check[6]) + 1)
                else:
                    if check[5] == "0":
                        check[2] = f[0]
                        check[3] = f[0]
                        check[4] = f[0]
                        check[5] = "1"
                        check[6] = str(int(check[6]) + 1)
                    else:
                        check[2] = f[0]
                        check[4] = f[0]
                        check[5] = str(int(check[5]) + 1)
                        check[6] = str(int(check[6]) + 1)
            elif len(f) == readcol+1:
                b.write("\t".join(check) + "\n")
                if f[1] == "+":
                    check = [chromosome,f[0],f[0],"-1","-1","0","1",f[readcol]]
                else:
                    check = [chromosome,f[0],f[0],f[0],f[0],"1","1",f[readcol]]
        c+=1
        print(c)
        d = a.readline()
    print(check)
    if int(check[3])>0:
        b.write("\t".join(check) + "\n")
    a.close()
    b.close()

def combinereadsinfolderx(ngsfold):
#10
    folderx = ngsfold + "bismarkfiles/bismarkmethyloutput/chroutput/CTOB/"
    finalfile = ngsfold + "bismarkfiles/bismarkmethyloutput/chroutput/CTOBreads.txt"
    a = os.listdir(folderx)
    c = 0
    while c < len(a):
        print(a[c])
        combinereadsinfilex(folderx + a[c], finalfile,3,a[c][:-4])
        c += 1
    folderx = ngsfold + "bismarkfiles/bismarkmethyloutput/chroutput/CTOT/"
    finalfile = ngsfold + "bismarkfiles/bismarkmethyloutput/chroutput/CTOTreads.txt"
    a = os.listdir(folderx)
    c = 0
    while c < len(a):
        print(a[c])
        combinereadsinfilex(folderx + a[c], finalfile,3,a[c][:-4])
        c += 1

def conversionsbyindcsinallrloops(inputfile,convcol,indccol, finalfile):
#11A cols are #of Cs, rows are # of conversions
    os.system("sort -k"+str(convcol+1)+","+str(convcol+1)+"n -k"+str(indccol+1)+","+str(indccol+1)+"n -f " + inputfile + " -T /home/jason/ -o " + inputfile)
    a = open(inputfile, "r")
    fl = open(finalfile, "w")
    b = a.readline()
    ic = 0
    c = 0
    currentconv = -1
    clist = []
    flist = []
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            if f[convcol] == currentconv:
                if clist.count(str(f[indccol])) == 1:
                    flist[clist.index(str(f[indccol]))] = str(int(flist[clist.index(str(f[indccol]))]) + 1)
                else:
                    clist.append(str(f[indccol]))
                    flist.append("1")
            else:
                print("b")
                if c == 0:
                    c += 1
                else:
                    fl.write(str(currentconv) + "\t" + "\t".join(flist) + "\n")
                count = 0
                flist = []
                currentconv = f[convcol]
                while count < len(clist):
                    flist.append("0")
                    count += 1
                print(flist)
                if clist.count(str(f[indccol])) == 1:
                    flist[clist.index(str(f[indccol]))] = str(int(flist[clist.index(str(f[indccol]))]) + 1)
                else:
                    clist.append(str(f[indccol]))
                    flist.append("1")
                print(flist)
                print(clist)
        b = a.readline()
    fl.write(str(currentconv) + "\t" + "\t".join(flist) + "\n")
    a.close()
    fl.close()
    a = open(finalfile, "r")
    b = a.readline()
    newlist = []
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            while len(f) < len(clist)+1:
                f.append("0")
            newlist.append(f)
        b = a.readline()
    a.close()
    fl = open(finalfile, "w")
    fl.write("-1" + "\t" + "\t".join(clist))
    c = 0
    while c < len(newlist):
        fl.write("\n" + "\t".join(newlist[c]))
        c += 1
    os.system("sort -k1,1n -f " + finalfile + " -T /home/jason/ -o " + finalfile)
    os.system("sort -k1,1n -f " + finalfile + " -T /home/jason/ -o " + finalfile)

def calculateavgconversion(datatable, maxconv):
#11Bi
    c = 1
    totalcs = 0
    convcs = 0
    while c < len(datatable):
        ci = 0
        conversionrow = datatable[c]
        while ci < len(conversionrow):
            if int(conversionrow[0])/int(datatable[0][ci]) < maxconv:
                convcs += int(conversionrow[ci])*int(conversionrow[0])
                totalcs += int(conversionrow[ci])*int(datatable[0][ci])
            ci += 1
        c += 1
    return float(convcs/totalcs)

def calculaterandomconversion(datatable):
#11B
    a = calculateavgconversion(datatable, 1.5)
    print(a)
    print(calculateavgconversion(datatable, a*2.5))
    return calculateavgconversion(datatable, a*2.5)

def inversefile(inputfile, outputfile):
#11C
    a = open(inputfile, "r")
    b = a.readline()
    final = []
    ic = 0
    e = b.replace("\n","")
    f = e.split("\t")
    if len(f)<2:
        ic += 1
    else:
        c = 0
        while c < len(f):
            final.append([f[c]])
            c += 1
    b = a.readline()
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            c = 0
            while c < len(f):
                final[c].append(f[c])
                c += 1
        b = a.readline()
    c = 0
    ff = open(outputfile, "w")
    while c < len(final):
        if c == 0:
            ff.write("\t".join(final[c]))
        else:
            ff.write("\n" + "\t".join(final[c]))
        c += 1
    ff.close()
    a.close()

def getinputdata(tstextfile):
#11Di
    tsfile = open(tstextfile, 'r')
    ts = tsfile.readlines()
    tsfile2 = [ter.replace('\n','') for ter in ts]
    tsfile3 = [ter.replace('"','') for ter in tsfile2]
    stsfile = [tsf.split('\t') for tsf in tsfile3]
    tsfile.close()
    return stsfile

def choose(x,n):
#11Dii1a
    final = 1
    if x > 0 and x < n:
        top = n
        c = n - 1
        while c > 0:
            top = top*c
            c -= 1
        bottom = x
        c = x - 1
        while c > 0:
            bottom = bottom*c
            c -= 1
        subt = n - x
        c = subt - 1
        while c > 0:
            subt = subt*c
            c -= 1
        final = top/(subt*bottom)
    return final

def binomial(p,x,n):
#11Dii1
    return choose(x,n) * p**x * (1 - p)**(n - x)


def expectedtablebinomial(probbychance,maxcs, finalf):
#11Dii
    c = 1
    ftable = [[]]
    intr = -1
    ff = open(finalf, "w")
    while intr < int(maxcs)+1:
        ftable[0].append(str(intr))
        intr += 1
    ff.write("\t".join(ftable[0]))
    while c < int(maxcs)+1:
        ftable.append([str(c)])
        n = 0
        while n < int(maxcs)+1:
            if n <= c:
                ftable[c].append(str(binomial(probbychance,n,c)))
            else:
                ftable[c].append(str(0))
            n += 1
        ff.write("\n" + "\t".join(ftable[c]))
        c += 1
    ff.close()

def rowreadscore(rowlist, rowprob):
#11Diii
    rowsum = 0
    ci = 1
    while ci < len(rowlist):
        rowsum += int(rowlist[ci])
        ci += 1
    ci = 2
    probsum = 1 - float(rowprob[1])
    csum = rowsum-int(rowlist[1])
    fl = [rowlist[0],"0"]
    newprob = 0
    while ci < len(rowlist) and float(rowprob[ci])*rowsum > 5:
        cvfr = int(rowlist[ci])/rowsum
        if cvfr > 0:
            newprob = (cvfr - float(rowprob[ci]))/cvfr
        if newprob < 0:
            newprob = 0
        fl.append(str(newprob))
        csum -= int(rowlist[ci])
        probsum -= float(rowprob[ci])
        ci += 1
    while ci < len(rowlist):
        cvfr = csum/rowsum
        if cvfr > 0:
            newprob = (cvfr-probsum)/cvfr
        if newprob < 0:
            newprob = 0
        fl.append(str(newprob))
        csum -= int(rowlist[ci])
        probsum -= float(rowprob[ci])
        ci += 1
    return fl

def combinelists(rowsuma, rowsumb, lista, listb):
#11Div
    weighta = int(rowsuma)/(int(rowsuma)+int(rowsumb))
    weightb = int(rowsumb)/(int(rowsuma)+int(rowsumb))
    final = [listb[0]]
    c = 1
    while c < len(listb):
        final.append(str(float(lista[c])*weighta + float(listb[c])*weightb))
        c += 1
    return final

def readscore(inputfile,probbychance,intfile, scorefile):
#11D
    a = getinputdata(inputfile)
    expectedtablebinomial(probbychance,a[-1][0], intfile)
    b = getinputdata(intfile)
    ff = open(scorefile, "w")
    ff.write("\t".join(a[0]))
    c = 1
    finalline = []
    while c < len(a):
        rowlist = a[c]
        rowsum = 0
        ci = 1
        while ci < len(rowlist):
            rowsum += int(rowlist[ci])
            ci += 1
        if rowsum > 1000:
            finalline = rowreadscore(rowlist, b[int(a[c][0])])
            rowsuma = rowsum
        if c < len(b):
            listb = rowreadscore(rowlist, b[int(a[c][0])])
            finalline = combinelists(rowsuma, rowsum, finalline, listb)
            rowsuma += rowsum
        ff.write("\n" + "\t".join(finalline))
        c += 1

def addscorestoinitialfile(inputfile,convcol,indccol, scorefile, finalfile):
#11E
    a = open(inputfile, "r")
    sc = getinputdata(scorefile)
    c = 0
    convrow = []
    while c < len(sc):
        convrow.append(sc[c][0])
        c += 1
    fl = open(finalfile, "w")
    b = a.readline()
    ic = 0
    c = 0
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            ai = convrow.index((f[indccol]))
            bi = sc[0].index(str(int(f[convcol])))
            f.append(sc[ai][bi])
            if c == 0:
                fl.write("\t".join(f))
            else:
                fl.write("\n" + "\t".join(f))
            c += 1
        b = a.readline()
    fl.close()
    a.close()

def addbisDRIPseqscorestofile(inputfile, finalfolder, stem):
#11
    print(inputfile)
    os.system("mkdir " + finalfolder)
    inta = finalfolder + stem + "inta.txt"
    intflipped = finalfolder + stem + "flipped.txt"
    intordered = finalfolder + stem + "flippedandordered.txt"
    intexpect = finalfolder + stem + "expected.txt"
    scorefile = finalfolder + stem + "score.txt"
    finalfile = finalfolder + stem + "withscores.txt"
    conversionsbyindcsinallrloops(inputfile,5,6, inta)
    a = getinputdata(inta)
    b = calculaterandomconversion(a)
    inversefile(inta,intflipped)
    os.system("sort -k1,1n " + intflipped + "> " + intordered)
    print("readscore start")
    readscore(intordered,b,intexpect,scorefile)
    print("addscorestoinitialfile start")
    addscorestoinitialfile(inputfile,5,6, scorefile, finalfile)
    print("finished runall")

def normreadscores(withscoresfile,normalizedfile):
#12
    a = open(withscoresfile, "r")
    ff = open(normalizedfile, "w")
    b = a.readline()
    ic = 0
    finalscore = 0
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            finalscore += float(f[-1])
        b = a.readline()
    a.close()
    a = open(withscoresfile, "r")
    b = a.readline()
    ic = 0
    c = 0
    while b != "" and ic < 100:
        e = b.replace("\n","")
        f = e.split("\t")
        if len(f)<2:
            ic += 1
        else:
            if c == 0 and f[-1]!="0" and float(f[-1])!=float("0"):
                ff.write("\t".join(f[:-1]) + "\t" + str((1000000*float(f[-1]))/finalscore))
                c += 1
            if c>0 and f[-1]!="0" and float(f[-1])!=float("0"):
                ff.write("\n" + "\t".join(f[:-1]) + "\t" + str((1000000*float(f[-1]))/finalscore))
        b = a.readline()
    a.close()
    ff.close()

def normscoresbyCs(filewnumreadsbytype, scorefile,finalfile):
#13Ai
    a = getinputdata(filewnumreadsbytype)
    s = getinputdata(scorefile)
    linescores = []
    linenums = []
    c = 1
    while c < len(a):
        cc = 1
        linescore = 0
        linenum = 0
        print(s[c])
        while cc < len(a[c]):
            linescore += int(a[c][cc])*float(s[c][cc])
            linenum += int(a[c][cc])
            cc += 1
        linescores.append(linescore)
        linenums.append(linenum)
        c += 1
    avgscore = sum(linescores)/sum(linenums)
    c = 1
    linenumssum = 0
    linescoressum = 0
    while c < len(s):
        if int(linenums[c-1])>1000:
            cc = 1
            multfactor = 1
            if linenums[c-1] > 0 and linescores[c-1] > 0:
                multfactor = avgscore*linenums[c-1]/linescores[c-1]
            while cc < len(a[c]):
                s[c][cc] = str(float(s[c][cc])*multfactor)
                cc += 1
            linenumssum = linenums[c - 1]
            linescoressum = linescores[c - 1]
            print(multfactor)
        else:
            linenumssum += linenums[c - 1]
            linescoressum += linescores[c - 1]
            cc = 1
            multfactor = 1
            if linenums[c - 1] > 0 and linescores[c - 1] > 0:
                multfactor = avgscore * linenumssum / linescoressum
            while cc < len(a[c]):
                s[c][cc] = str(float(s[c][cc]) * multfactor)
                cc += 1
        c += 1
    ff = open(finalfile,"w")
    ff.write("\t".join(s[0]))
    c = 1
    while c < len(s):
        ff.write("\n" + "\t".join(s[c]))
        c += 1
    ff.close()

def normreadscorebyCs(inputfile,probbychance,intfile, scorefile,normscorefile):
#13A
    a = getinputdata(inputfile)
    expectedtablebinomial(probbychance,a[-1][0], intfile)
    b = getinputdata(intfile)
    ff = open(scorefile, "w")
    ff.write("\t".join(a[0]))
    c = 1
    finalline = []
    while c < len(a):
        rowlist = a[c]
        rowsum = 0
        ci = 1
        while ci < len(rowlist):
            rowsum += int(rowlist[ci])
            ci += 1
        if rowsum > 1000:
            finalline = rowreadscore(rowlist, b[c])
            rowsuma = rowsum
        if c < len(b):
            listb = rowreadscore(rowlist, b[c])
            finalline = combinelists(rowsuma, rowsum, finalline, listb)
            rowsuma += rowsum
        ff.write("\n" + "\t".join(finalline))
        c += 1
    ff.close()
    normscoresbyCs(inputfile, scorefile, normscorefile)

def runallnormbyC(inputfile, finalfolder, stem):
#13
    print(inputfile)
    os.system("mkdir " + finalfolder)
    inta = finalfolder + stem + "inta.txt"
    intflipped = finalfolder + stem + "flipped.txt"
    intordered = finalfolder + stem + "flippedandordered.txt"
    intexpect = finalfolder + stem + "expected.txt"
    intscorefile = finalfolder + stem + "prescore.txt"
    scorefile = finalfolder + stem + "score.txt"
    finalfile = finalfolder + stem + "withscores.txt"
    normreadscorebyCs(intordered,b,intexpect,intscorefile,scorefile)
    print("addscorestoinitialfile start")
    addscorestoinitialfile(inputfile,5,6, scorefile, finalfile)
    print("finished runall")

def wbedline(filea,chro,loc,value1, value2):
#14Ai
    value3 = 0
    if value2 > 0:
        value3 = value1/value2
    filea.write("\n" + chro + "\t" + str(loc) + "\t" + str(value1) + "\t" + str(value3) + "\t" + str(value2))

def addtobed(templ,chro, file1):
#14A
    templ.sort()
    a = []
    if templ[0][1] == "-":
        a.append([templ[0][0],1,1])
    else:
        a.append([templ[0][0],0,1])
    countera = 1
    while countera < len(templ):
        if templ[countera][0] == a[-1][0]:
            if templ[countera][1] == "-":
                a[-1][1] += 1
                a[-1][2] += 1
            else:
                a[-1][2]+=1
            if countera == len(templ)-1:
                wbedline(file1,chro,a[-1][0],a[-1][1],a[-1][2])
                countera += 1
            else:
                countera += 1
        else:
            wbedline(file1,chro,a[-1][0],a[-1][1],a[-1][2])
            if templ[countera][1] == "-":
                a.append([templ[countera][0],1,1])
            else:
                a.append([templ[countera][0],0,1])
            if countera == len(templ)-1:
                wbedline(file1,chro,a[-1][0],a[-1][1],a[-1][2])
                countera += 1
            else:
                countera += 1

def makeconversionsfile(chrfolder, outputfile):
#14
    start = os.listdir(chrfolder)
    file1 = open(outputfile + "conversionsperc.txt", "w")
    counter = 0
    while counter < len(start):
        tf = chrfolder + start[counter]
        templ = getinputdata(tf)
        addtobed(templ, start[counter][:-4], file1)
        counter += 1

def runbisDRIPseqpipeline(ngsfold,outputfolder,name):
    runflexbar(ngsfold)
    runbismark(ngsfold)
    prepareunmappedfull(ngsfold)
    secondrunbismark(ngsfold)
    finalrunbismark(ngsfold)
    runmethylationextractor(ngsfold)
    rundepulication(ngsfold)
    actualextract(ngsfold)
    executedatadump(ngsfold)
    combinereadsinfolderx(ngsfold)
    inputfile = ngsfold + "bismarkfiles/bismarkmethyloutput/chroutput/CTOBreads.txt"
    addbisDRIPseqscorestofile(inputfile, outputfolder, name + "CTOB")
    normreadscores(outputfolder + name + "CTOBwithscores.txt",outputfolder + name + "CTOBstdnormalized.txt")
    runallnormbyC(inputfile, outputfolder, name + "CTOB")
    inputfile = ngsfold + "bismarkfiles/bismarkmethyloutput/chroutput/CTOTreads.txt"
    addbisDRIPseqscorestofile(inputfile, outputfolder, name + "CTOT")
    normreadscores(outputfolder + name + "CTOBwithscores.txt",outputfolder + name + "CTOTstdnormalized.txt")
    runallnormbyC(inputfile, outputfolder, name + "CTOT")
    makeconversionsfile(ngsfold + "bismarkfiles/bismarkmethyloutput/chroutput/CTOB/", outputfolder + name + "CTOBconversions.txt")
    makeconversionsfile(ngsfold + "bismarkfiles/bismarkmethyloutput/chroutput/CTOT/", outputfolder + name + "CTOTconversions.txt")

