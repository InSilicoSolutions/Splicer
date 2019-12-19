import sqlite3
import os
import argparse

parser = argparse.ArgumentParser(description='geneModelImport reads reference transcript models into the Splicer database.  It accepts either a GTF file or a set of  UCSC knownGene files.', usage='geneModelImport database_directory  GTF|UCSC GTF_File|knownGene.txt kgTxInfo.txt kgXRef.txt')
parser.add_argument("Database", help='The path to where you want to store the database file.')
parser.add_argument("Type", choices=['GTF', 'UCSC'], help='The type of the files to be imported.')
parser.add_argument("Input", nargs='+', help='GTF input file or UCSC input files.')
args = parser.parse_args()

    
#class that holds the data of an exon
#these are made for each exon line in a transcript, updated with the cds and utr info, and then inserted into the database
class Exon:
    def __init__(self, tid, eid, exonNum, chrom, startPos, stopPos):
        self.tid = tid
        self.eid = eid
        self.exonNum = exonNum
        self.chrom = chrom
        self.startPos = startPos
        self.stopPos = stopPos
        self.cdsStart = None
        self.cdsStop = None
        self.utrStart = None
        self.utrStop = None
    
    #inserts this exon into the database
    def insert(self):
        if self.utrStart != None:
            #if there is more than one utr join them with a dash
            if len(self.utrStart) == 2 :
                self.utrStart = "-".join(self.utrStart)
                self.utrStop = "-".join(self.utrStop)
            #if there is only one utr take it out of the list so it is just a string
            elif len(self.utrStart) == 1 :
                self.utrStart = self.utrStart[0]
                self.utrStop = self.utrStop[0]
        c.execute("INSERT INTO Exon VALUES("+str(self.tid)+", "+str(self.eid)+", "+str(self.exonNum)+", '"+str(self.chrom)+"', "+str(self.startPos)+", "+str(self.stopPos)+", '"+str(self.cdsStart)+"', '"+str(self.cdsStop)+"', '"+str(self.utrStart)+"', '"+str(self.utrStop)+"')")

    #checks if a given utr start is within the bounds of this exon
    def contains(self, utrStart):
        if utrStart >= self.startPos and utrStart <= self.stopPos:
            return True
        else:
            return False

#get connection to the sqlite database
conn = sqlite3.connect(args.Database + os.path.sep + 'splice.sqlite', isolation_level=None)
c = conn.cursor()

#create or rebuild tables
#gene table (Index, Gene Symbol, Chromosome, Strand)
c.execute("DROP TABLE IF EXISTS Gene;")
c.execute('''CREATE TABLE Gene
             (Gene_ID INTEGER PRIMARY KEY NOT NULL DEFAULT NULL,
              Reference_ID INTEGER NOT NULL DEFAULT NULL,
              Symbol varchar(20) NOT NULL DEFAULT NULL, 
              Chromosome varchar(30) NOT NULL DEFAULT NULL,
              Strand varchar(1) NOT NULL DEFAULT NULL,
              Start_Position INTEGER NOT NULL DEFAULT NULL,
              Stop_Position INTEGER NOT NULL DEFAULT NULL,
              Type varchar(30) DEFAULT NULL,
              Source varchar(30) DEFAULT NULL);''')
c.execute("CREATE INDEX idx_Gene_Symbol ON Gene(symbol);")

c.execute("DROP TABLE IF EXISTS Transcript;")
c.execute('''CREATE TABLE Transcript
             (Gene_ID INTEGER NOT NULL DEFAULT NULL,
              Transcript_ID INTEGER PRIMARY KEY NOT NULL DEFAULT NULL,
              Transcript_Reference_ID INTEGER NOT NULL DEFAULT NULL,
              Chromosome varchar(30) NOT NULL DEFAULT NULL,
              Start_Position INTEGER NOT NULL DEFAULT NULL,
              Stop_Position INTEGER NOT NULL DEFAULT NULL,
              Type varchar(30) NOT NULL DEFAULT NULL);''')
c.execute("CREATE INDEX idx_Transcript_Gene ON Transcript(Gene_ID);")

c.execute("DROP TABLE IF EXISTS Exon;")
c.execute('''CREATE TABLE Exon
             (Transcript_ID INTEGER NOT NULL DEFAULT NULL,
              Exon_ID INTEGER PRIMARY KEY NOT NULL DEFAULT NULL,
              Exon_Number INTEGER NOT NULL DEFAULT NULL,
              Chromosome varchar(30) NOT NULL DEFAULT NULL,
              Start_Position INTEGER NOT NULL DEFAULT NULL,
              Stop_Position INTEGER NOT NULL DEFAULT NULL,
              CDS_Start varchar(30) DEFAULT NULL,
              CDS_Stop varchar(30) DEFAULT NULL,
              UTR_Start varchar(30) DEFAULT NULL,
              UTR_Stop varchar(30) DEFAULT NULL);''')
c.execute("CREATE INDEX idx_Exon_Transcript ON Exon(Transcript_ID);")


#column 8 of gtf format is s dictionary in string form
#this function makes that string into an actual dictionary
def dictGen(stringDict):
    retDict = {}
    sdSplitList = stringDict.split(";")
    #pop off the empty string at the end of the list
    sdSplitList.pop()
    #loop through remaining strings and add them to the dict
    for entry in sdSplitList:
        eSplit = entry.replace('"','').split(" ")
        if eSplit[0] != "":
            retDict[eSplit[0]] = eSplit[1]
        else:
            retDict[eSplit[1]] = eSplit[2]
    return retDict

#evaluates the CDS and UTR for a coding ucsc exon
def UcscCoding(started, ended, codingTrans, lList, exonStart, exonStop, cdsStart, cdsStop):
    exonCdsStart = None
    exonCdsStop = None
    exonUtrStart = None
    exonUtrStop = None
    if(codingTrans[lList[0]]):
        #whole exon is utr
        if not started or ended:
            exonUtrStart = exonStart
            exonUtrStop = exonStop
                
        #cds start is in the exon
        if cdsStart <= exonStop and cdsStart >= exonStart:
            started = True
            if cdsStart == exonStart:
                exonCdsStart = exonStart
            else:
                exonUtrStart = str(exonStart)
                exonUtrStop = str(int(cdsStart)-1)
                exonCdsStart = cdsStart
                        
        #only look for the end if it is started
        if started:
            #if the cds did not start this exon then the cds start is the exon start
            if exonCdsStart == None:
                exonCdsStart = exonStart
            #if cds stop is in the exon
            if cdsStop <= exonStop and cdsStop >= exonStart:
                ended = True
                if cdsStop == exonStop:
                    exonCdsStop = exonStop
                else:
                    exonCdsStop = cdsStop
                    #if there is already a utr add a second one after a dash
                    if exonUtrStart == None:
                        exonUtrStart = str(int(cdsStop)+1)
                        exonUtrStop = str(exonStop)
                    else:
                        exonUtrStart = str(exonUtrStart) + "-" + str(int(cdsStop)+1)
                        exonUtrStop = str(exonUtrStop) + "-" + str(exonStop)
            else:
                exonCdsStop = exonStop
                
    return   {"cdsStart":exonCdsStart, "cdsStop":exonCdsStop, "utrStart":exonUtrStart, "utrStop":exonUtrStop, "started":started, "ended":ended}      

#populate tables from a gtf input file
def processGTF():
    #open the input file    INPUT FILE FOR TESTING: Homo_sapiens.GRCh37.87.gtf
    infile = open(args.Input[0], 'r')    
    geneIndex = 0
    tranIndex = 0
    exonIndex = 0
    exonList = []
    elTracker = 0
    col8 = {}
    linecount = 0
    c.execute("begin")
    for line in infile:
        linecount+=1
        print(linecount)
        lList = line.replace('\n','').split("\t")
        if len(lList)>2:
                    
            if lList[2] == "gene":
                geneIndex += 1
                col8 = dictGen(lList[8])
                c.execute("INSERT INTO Gene VALUES("+str(geneIndex)+", '"+col8['gene_id']+"', '"+col8['gene_name']+"', '"+lList[0]+"', '"+lList[6]+"', '"+lList[3]+"', '"+lList[4]+"', '"+col8['gene_biotype']+"', '"+col8['gene_source']+"')")
        
            elif lList[2] == "transcript":
                #add exon stuff to the database and reset it
                for exon in exonList:
                    exon.insert()
                elTracker = 0
                exonList.clear()
                #insert transcript stuff    
                tranIndex += 1
                col8 = dictGen(lList[8])
                c.execute("INSERT INTO Transcript VALUES("+str(geneIndex)+", "+str(tranIndex)+", '"+col8['transcript_id']+"', '"+lList[0]+"', '"+lList[3]+"', '"+lList[4]+"', '"+col8['transcript_biotype']+"')")
            #if this line is a new exon create a new exon object
            elif lList[2] == "exon":
                exonIndex += 1
                col8 = dictGen(lList[8])
                exonList.append(Exon(str(tranIndex), str(exonIndex), col8['exon_number'], lList[0], lList[3], lList[4]))
                
            elif lList[2] == "CDS":
                exonList[-1].cdsStart = lList[3]
                exonList[-1].cdsStop = lList[4]
                
            elif lList[2].find("_prime_utr") != -1:
                for i in range(elTracker, len(exonList)):
                    if exonList[i].contains(lList[3]):
                        elTracker = i
                        if exonList[i].utrStart == None:
                            exonList[i].utrStart = [lList[3]]
                            exonList[i].utrStop = [lList[4]]
                            break
                        else:
                            exonList[i].utrStart.append(lList[3])
                            exonList[i].utrStop.append(lList[4])
                            break
    
    #add the exons from the last transcript in the file
    if len(exonList) != 0:
        for exon in exonList:
            exon.insert()
    
    c.execute("commit")
    print("done")

#populate tables from ucsc input files  
def processUCSC():
    #read kgXref and build a dictionary(enst -> gene symbol)
    infile = open(args.Input[2], 'r')    
    enstSymbol = {}
    for line in infile:
        lList = line.replace('\n','').split("\t")
        enstSymbol[lList[0]] = lList[4]
    #read kgTxInfo to make and dictionary about which transcripts are coding
    infile = open(args.Input[1], 'r')    
    codingTrans = {}
    typeTrans = {}
    for line in infile:
        lList = line.replace('\n','').split("\t")
        typeTrans[lList[0]] = lList[1]
        if lList[1] == 'coding' and lList[11] == '1' and lList[12] == '1':
            codingTrans[lList[0]] = True
        else:
            codingTrans[lList[0]] = False
    
    #populate gene table
    #open input file       TESTING FILE: knownGene.txt
    infile = open(args.Input[0], 'r')
    geneIdx = 1
    tranIdx = 0
    exonIdx = 1
    currentSymbol = ""
    currentChrom = ""
    currentStrand = ""
    geneCoding = False
    maxStop = None
    minStart = None
    c.execute("begin")
    for line in infile:
        lList = line.replace('\n','').split("\t")
        starts = lList[8].split(",")
        stops = lList[9].split(",")
        #add one to all starts and cds start
        for i in range(len(starts)-1):
            #add one to start
            starts[i] = str(int(starts[i]) + 1)
        #add one to cds start
        lList[5] = str(int(lList[5]) + 1)
        #add one to gene start
        lList[3] = str(int(lList[3]) + 1)
        
        symbol = enstSymbol[lList[0]]
        if currentSymbol=="":
            currentSymbol = symbol 
            currentChrom = lList[1]
            currentStrand = lList[2]
            minStart = lList[3]
            maxStop = lList[4]
        #check if the line belongs to the current symbol
        if symbol != currentSymbol :
            
            geneType = "noncoding"
            if geneCoding:
                geneType = "coding"
            c.execute("INSERT INTO gene VALUES("+str(geneIdx)+", '""', '"+currentSymbol+"', '"+currentChrom+"', '"+currentStrand+"', '"+minStart+"', '"+maxStop+"', '"+geneType+"', 'UCSC KnownGene')")
            geneCoding = False
            currentSymbol = symbol 
            currentChrom = lList[1]
            currentStrand = lList[2]
            minStart = lList[3]
            maxStop = lList[4]
            print(geneIdx)
            geneIdx+=1
        #always insert the transcript info
        tranIdx+=1
        c.execute("INSERT INTO Transcript VALUES("+str(geneIdx)+", '"+str(tranIdx)+"', '"+lList[0]+"', '"+lList[1]+"', "+lList[3]+", "+lList[4]+", '"+typeTrans[lList[0]]+"')")
        #if one transcript is coding set type to coding for this gene
        if codingTrans[lList[0]]:
            geneCoding = True
        
        minStart = minStart if lList[3] > minStart else lList[3]
        maxStop = maxStop if lList[4] < maxStop else lList[4]
        started = False
        ended = False
        for i in range(len(starts)-1):
            #find out exon number (counts backwards if negative strand)
            iNumber = i+1
            exonNumber = iNumber if lList[2] == '+' else len(starts)-iNumber
            #if its coding figure out CDS and UTR
            retCoding = UcscCoding(started, ended, codingTrans, lList, starts[i], stops[i], lList[5], lList[6])
            #update started and ended trackers
            started = retCoding['started']
            ended = retCoding['ended']
            c.execute("INSERT INTO Exon VALUES("+str(tranIdx)+", "+str(exonIdx)+", "+str(exonNumber)+", '"+lList[1]+"', "+starts[i]+", "+stops[i]+", '"+str(retCoding['cdsStart'])+"', '"+str(retCoding['cdsStop'])+"', '"+str(retCoding['utrStart'])+"', '"+str(retCoding['utrStop'])+"')")
            exonIdx+=1
    #insert the last gene
    geneType = "noncoding"
    if geneCoding:
        geneType = "coding"
    c.execute("INSERT INTO gene VALUES("+str(geneIdx)+", '""', '"+currentSymbol+"', '"+currentChrom+"', '"+currentStrand+"', '"+minStart+"', '"+maxStop+"', '"+geneType+"', 'UCSC KnownGene')")
    print(geneIdx)        
    c.execute("commit")
    print('done')
    
if args.Type == "GTF":
    processGTF()
elif args.Type == "UCSC":
    if (len(args.Input) < 3):
        parser.print_help()
    else:
        processUCSC()
    
        
