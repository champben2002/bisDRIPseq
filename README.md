# bisDRIPseq
This project includes files required to perform bisDRIP-seq analysis

(1) processingbisDRIPseqreads.py generates bisDRIP-seq scores for each bisDRIP-seq read from a fastq file.

(2) bisDRIPseqmetaplotanalysis.py generates metaplots of bisDRIP-seq scores relative to a set of features.

(3) regionbisDRIPseqscores.py calculates the bisDRIP-seq score of regions in the genome from bisDRIP-seq read scores.

(4) Monte_Carlo_random_assign_reads_to_regions.py performs a Monte Carlo simulation in which reads are assigned to regions at random.

(5) Monte_Carlo_for_shuffling_bisDRIPseq_scores.py randomly shuffles bisDRIP-seq scores to reads.


(1) processingbisDRIPseqreads.py

  (A) Dependencies:
  
    (i) flexbar must be installed
    
      (a) Note: change '/opt/flexbar_v2.5_linux64/flexbar' on line 86 to the location of your flexbar installation
      
      (b) Note: change '/opt/flexbar_v2.5_linux64/flexbar' on line 162 to the location of your flexbar installation
      
    (ii) bismark must be installed as well as its dependencies
    
      (a) Note: change '/home/jason/bismark_v0.14.3/bismark' on line 112 to the location of your bismark installation
      
      (b) Note: change '/home/jason/bismark_v0.14.3/bismark' on line 142 to the location of your bismark installation
      
      (c) Note: change '/home/jason/bismark_v0.14.3/bismark' on line 179 to the location of your bismark installation
      
      (d) Note: change '/home/jason/bismark_v0.14.3/bismark_methylation_extractor' on line 224 to the location of your bismark   installation
      
      (e) Note: change '/home/jason/bismark_v0.14.3/bismark_methylation_extractor' on line 225 to the location of your bismark installation
      
      (f) Note: change '/home/jason/bismark_v0.14.3/deduplicate_bismark' on line 233 to the location of your bismark installation
      
      (g) Note: change '/home/jason/bismark_v0.14.3/deduplicate_bismark' on line 240 to the location of your bismark installation
      
      (h) Note: change '/home/jason/bismark_v0.14.3/bismark_methylation_extractor' on line 268 to the location of your bismark installation
      
      (i) Note: change '/home/jason/bismark_v0.14.3/bismark_methylation_extractor' on line 269 to the location of your bismark installation
      
      (g) Note: change '/home/jason/bismark_v0.14.3/deduplicate_bismark' on line 233 to the location of your bismark installation
      
      (h) Note: change '/home/jason/bismark_v0.14.3/deduplicate_bismark' on line 233 to the location of your bismark installation
  
  
  (B) Run program by running the algorithm runbisDRIPseqpipeline()
  
  
  (C) Output of runbisDRIPseqpipeline(): The output is a set of files in the directory OUTPUTFOLDER as defined by the user
    
    (i) Standard normalized bisDRIP-seq read files written to the directory "outputfolder" with the file names NAME + "CTOBstdnormalized.txt" and NAME + "CTOTstdnormalized.txt" where NAME is a user inputed variable. The contents of these file will contain bisDRIP-seq reads with their associated bisDRIP-seq scores. Each row of the files will contain information regarding a single read. The first column of each row will contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row will contain the start position of each read. The third column in each row will contain the end position of each read. The final column of each row will contain the bisDRIP-seq score. Each of these two files represents a different strand.
    
    (ii) Cytosine normalized bisDRIP-seq read files written to the directory "outputfolder" with the file names NAME + "CTOBcnormwithscores.txt" and NAME + "CTOTcnormwithscores.txt" where NAME is a user inputed variable. The contents of these file will contain bisDRIP-seq reads with their associated bisDRIP-seq scores. Each row of the files will contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row will contain the start position of each read. The third column in each row will contain the end position of each read. The final column of each row will contain the bisDRIP-seq score. Each of these two files represents a different strand.
    
    (iii) Raw bisDRIP-seq read files written to the directory "outputfolder" with the file names NAME + "CTOBwithscores.txt" and NAME + "CTOTwithscores.txt" where NAME is a user inputed variable. The contents of these file will contain bisDRIP-seq reads with their associated bisDRIP-seq scores. Each row of the files will contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row will contain the start position of each read. The third column in each row will contain the end position of each read. The final column of each row will contain the bisDRIP-seq score. Each of these two files represents a different strand.
    
    (iv) BisDRIP-seq conversion files written to the directory "outputfolder" with the file names NAME + "CTOBconversions.txt" and NAME + "CTOTconversions.txt" where NAME is a user inputed variable. The contents of these file will contain each cytosine nucleotide position in the genome with the number of reads and conversions associated with that cytosine. Each row of the files will contain information regarding a single cytosine position in the genome. The first column of each row will contain the chromosome of the cytosine using the notation chrA where A is the number or letter of the chromosome. The second column of each row will contain the position of each cytosine. The third column in each row will contain the number of reads where the position was read as a T (ie converted). The fourth column of each row will contain the fraction of reads where the position was read as a T (ie the fraction of times the cytosine was converted). The final column will  of each row will contain the total number of reads that aligned to the position. Each of these two files represents a different strand.
    
    (v) Other files: Other folders and files contain the various intermediary steps to produce these output files. These files are not deleted by default.
 
 
  (C) Input of runbisDRIPseqpipeline():
  
     (i) ngsfold: a fastq file of bisDRIP-seq reads generated by CASAVA 1.8.2
     
     (ii) outputfolder: directory where output will be written
     
     (iii)name: stem of the name for output files


(2) bisDRIPseqmetaplotanalysis.py

  (A) Run program by running the algorithm allproximalscorestosinglentfeatures()


  (B) Output: The output is a multi-column tab-delimited file. In the first column are the positions from -disfromfeature to +disfromfeature, with zero indicating the location of the feature(s). In each subsequent column is the sum of the bisDRIP-seq scores at that position across all features for a given bisDRIP-seq sample.


  (C) Input of allproximalscorestosinglentfeatures():

    (i) featurefile: a tab-delimited file that contains the location of each feature. The first row must be a header. Each following row must contain a feature and each feature must be on a seperate row. The first column of each feature row must contain the chromosome of the feature using the notation chrA where A is the number or letter of the chromosome. The second column of each feature row must contain the genomic position of the feature.
    
    (ii) inputfolder: a directory containing a set of bisDRIP-seqs read score files. Each of these files should be for a different bisDRIP-seq sample. The contents of the file should contain bisDRIP-seq reads with their associated bisDRIP-seq scores. Each row of the files should contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row should contain the start position of each read. The third column in each row should contain the end position of each read. The final column of each row should contain the bisDRIP-seq score. 
    
    (iii) outputfolder: directory where output should be written
    
    (iv) name: the output will be written to outputfolder/ + name + final.txt
    
    (v) direction: if "negative", then the metaplot will be flipped around so that nucleotides upstream of the feature will be viewed as +X bp from the feature. This is valuable for genes which can be orientated in either direction. Must otherwise be "positive".
    
    (vi) disfromfeature: This determines how far from the feature the metaplot will extend from features. Note that disfromfeature must be smaller than the minimum distance between any two features


(3) regionbisDRIPseqscores.py 
  
  (A) Run program by running the algorithm proximalscorestart()
  
  
  (B) Output of proximalscorestart(): The output is a multi-column tab-delimited file. Each row contains a different feature of interest. In the first column of each row is the chromosome of the feature. In the second column of each row is the position of the feature. In the third column of each row is the bisDRIP-seq score associated with that feature.
  
  
  (C) Input of proximalscorestart():
  
    (i) featurefile: a tab-delimited file that contains the location of each feature. The first row must be a header. Each following row must contain a feature and each feature must be on a seperate row. The first column of each feature row must contain the chromosome of the feature using the notation chrA where A is the number or letter of the chromosome. The second column of each feature row must contain the genomic position of the feature.
    
    (ii) inputfile: a bisDRIP-seqs read score file. The contents of the file should contain bisDRIP-seq reads with their associated bisDRIP-seq scores. Each row of the files should contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row should contain the start position of each read. The third column in each row should contain the end position of each read. The final column of each row should contain the bisDRIP-seq score. 
    
    (iii) folder - directory where output should be written
    
    (iv) name - output will be written into the file: folder + name + "final.txt"
    
    (v) Direction: if "positive" than the output provides the score after and including the feature otherwise you get back the score before and including the feature
    
    (vi) disfromfeature: The returned score will be of the sum of the bisDRIP-seq scores in the region between the feature and the disfromfeature


(4) Monte_Carlo_random_assign_reads_to_regions.py
  
  (A) Run program by running the algorithm runrandomizeronafileRRF
  
  
  (B) Output of runrandomizeronafileRRF(): A multi-column tab-delimited file with bisDRIP-seq reads from a single chromosome associated with bisDRIP-seq scores and a randomly assigned location. All reads in this file are associated with scores above zero. Each row of the file will contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row should contain the randomly-assigned start position of each read. The third column in each row should contain the assigned end position of each read. The final column of each row should contain the bisDRIP-seq score of the read. 

  
  (C) Input of runrandomizeronafileRRF():
    
    (i) readfile:  A multi-column tab-delimited file with bisDRIP-seq reads from a single chromosome associated with bisDRIP-seq scores. Reads associated with bisDRIP-seq scores that equal zero may be excluded from this file. Each row of the file should contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row should contain the start position of each read. The third column in each row should contain the end nucleotide of each read. The final column of each row should contain the bisDRIP-seq score of the read.
    
    (ii) checkfile:  A multi-column tab-delimited file with regions from a single chromosome that contain bisDRIP-seq reads. Each row of the file should contain information regarding a single region. The first column of each row must contain the chromosome of the region using the notation chrA where A is the number or letter of the chromosome. The second column of each region should contain the median positon of the region.
    
    (iii) tempstem: a directory and possible start of file name where files will be temporary written and read from.
    
    (iv) outputfile: name (including directory) of the file that the contents will be written into. 
    
    
(5) Monte_Carlo_for_shuffling_bisDRIPseq_scores.py
  
  (A) Run program by running the algorithm runswitchlabelsonchrms()
  
  
  (B) Output of runswitchlabelsonchrms(): A multi-column tab-delimited file with bisDRIP-seq reads from a single chromosome associated with randomly assigned bisDRIP-seq scores. All reads in this file are associated with scores above zero. Each row of the file will contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row should contain the start position of each read. The third column in each row should contain the end position of each read. The final column of each row should contain the randomly-assigned bisDRIP-seq score. 
  
  
  (C) Input of runswitchlabelsonchrms():
    
    (i) readbisDRIPseqscorefile:  A multi-column tab-delimited file with bisDRIP-seq reads from a single chromosome associated with bisDRIP-seq scores. Reads associated with bisDRIP-seq scores that equal zero may be excluded from this file. Each row of the file should contain information regarding a single read. The first column of each row must contain the chromosome of the read using the notation chrA where A is the number or letter of the chromosome. The second column of each row should contain the start position of each read. The third column in each row should contain the end position of each read. The final column of each row should contain the bisDRIP-seq score of the read.
    
    (ii) allreadfile:  A multi-column tab-delimited file with all bisDRIP-seq reads from a single chromosome. Each row of the file should contain information regarding a single read. Each row should have multiple columns.
    
    (iii) outputfolder: a directory where the output will be written into
    
    (iv) name - name of the file that the contents will be written into. 
