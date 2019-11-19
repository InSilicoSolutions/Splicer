import sys
import sqlite3
import os
from _sqlite3 import IntegrityError

class Exon:
    def __init__(self, tid, eid, exonNum, chrom, startPos, stopPos):
        self.tid = tid
        self.eid = eid
        self.exonNum = exonNum
        self.chrom = chrom
        self.startPos = startPos
        self.stopPos = stopPos
        self.cdsStart = 0
        self.cdsStop = 0
        self.utrStart = None
        self.utrStop = None
    
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
        c.execute("INSERT INTO Exon VALUES("+str(self.tid)+", "+str(self.eid)+", "+str(self.exonNum)+", '"+str(self.chrom)+"', "+str(self.startPos)+", "+str(self.stopPos)+", "+str(self.cdsStart)+", "+str(self.cdsStop)+", '"+str(self.utrStart)+"', '"+str(self.utrStop)+"')")

        
    def contains(self, utrStart):
        if utrStart >= self.startPos and utrStart <= self.stopPos:
            return True
        else:
            return False
"""        
class SubExon:
    def __init__(self, name, start, stop):
        self.name = name
        self.start = start
        self.stop = stop
        
        
class exonDataContainer:
    def __init__(self, starts, stops, cdsStart, cdsStop, strand):
        self.starts = starts
        self.stops = stops
        self.cdsStarts = [cdsStart]
        self.cdsStops = [cdsStop]
        self.strand = strand
    
    #add a new line of data to the variables
    def addLine(self, starts, stops, cdsStart, cdsStop):
        self.starts.update(starts)
        self.stops.update(stops)
        self.cdsStarts.append(cdsStart)
        self.cdsStops.append(cdsStop)
        
    def datInsert(self):
        #add all the positions to a list and sort them
        posList = list(self.starts)
        posList.extend(self.stops)
        posList.sort()
        #if negative strand reverse
        if self.strand == "-":
            temp = self.starts
            self.starts = self.stops
            self.stops = temp
            posList.reverse()
        #boolean used to keep track of if the loop has seen a stop position
        boolStop = False
        current = SubExon("1.1", posList[0], None)
        exonCount = 1
        subExonCount = 1
        exonList = []
        #loop through the positions making appropriate subExon objects
        for i in range(1,len(posList)):
            
            #note the positions here are splice starts and stops so where a splice stops is where an exon starts
            if posList[i] in self.stops or boolStop == False:
                #set the stop marker if this is a stop position
                if posList[i] in self.stops:
                    boolStop = True
                #update the current sub exon with its end position
                current.stop = posList[i]
                #if this is the first time the exon has been split change the name of the exon to have a .1
                subExonCount += 1
                #add the previous exon to the list before writing over the current variable
                exonList.append(current)
                #make the new sub exon object
                current = SubExon(str(exonCount)+"."+str(subExonCount), posList[i], None)
                
            else:
                boolStop = False
                #if the exon has no subexons remove the '.1' from the name
                if current.name.split('.')[1] == "2":
                    exonList[-1].name = exonList[-1].name.split(".")[0]
                exonCount += 1
                subExonCount = 1
                current = SubExon(str(exonCount)+"."+str(subExonCount), posList[i], None)
        if current.name.split('.')[1] == "2":
            exonList[-1].name = exonList[-1].name.split(".")[0]
            #for exon in exonList:
            #add to database   
        print(posList)

"""
#get connection to the sqlite database
conn = sqlite3.connect("E:\speedSplice" + os.path.sep + 'splice.sqlite', isolation_level=None)
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
c.execute("CREATE INDEX gSymbol ON Gene(symbol);")

c.execute("DROP TABLE IF EXISTS Transcript;")
c.execute('''CREATE TABLE Transcript
             (Gene_ID INTEGER NOT NULL DEFAULT NULL,
              Transcript_ID INTEGER PRIMARY KEY NOT NULL DEFAULT NULL,
              Transcript_Reference_ID INTEGER NOT NULL DEFAULT NULL,
              Chromosome varchar(30) NOT NULL DEFAULT NULL,
              Start_Position INTEGER NOT NULL DEFAULT NULL,
              Stop_Position INTEGER NOT NULL DEFAULT NULL);''')

c.execute("DROP TABLE IF EXISTS Exon;")
c.execute('''CREATE TABLE Exon
             (Transcript_ID INTEGER NOT NULL DEFAULT NULL,
              Exon_ID INTEGER PRIMARY KEY NOT NULL DEFAULT NULL,
              Exon_Number INTEGER NOT NULL DEFAULT NULL,
              Chromosome varchar(30) NOT NULL DEFAULT NULL,
              Start_Position INTEGER NOT NULL DEFAULT NULL,
              Stop_Position INTEGER NOT NULL DEFAULT NULL,
              CDS_Start INTEGER DEFAULT NULL,
              CDS_Stop INTEGER DEFAULT NULL,
              UTR_Start varchar(30) DEFAULT NULL,
              UTR_Stop varchar(30) DEFAULT NULL);''')

#turn a string of a dictionary into a dictionary
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

if sys.argv[0] == "gtf":
    infile = open('Homo_sapiens.GRCh37.87.gtf', 'r')    
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
                c.execute("INSERT INTO Transcript VALUES("+str(geneIndex)+", "+str(tranIndex)+", '"+col8['transcript_id']+"', '"+lList[0]+"', '"+lList[3]+"', '"+lList[4]+"')")
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
elif sys.argv[0] == "ucsc":
    #read kgXref and build a dictionary(enst -> gene symbol)
    infile = open('kgXref.txt', 'r')    
    enstSymbol = {}
    for line in infile:
        lList = line.replace('\n','').split("\t")
        enstSymbol[lList[0]] = lList[4]
    
    #populate gene table
    infile = open('knownGene.txt', 'r')
    geneIdx = 1;
    tranIdx = 1;
    currentSymbol = ""
    currentChrom = ""
    currentStrand = ""
    maxStop = None
    minStart = None
    exonData = {} #symbol -> exon data object
    c.execute("begin")
    for line in infile:
        lList = line.replace('\n','').split("\t")
        starts = lList[8].split(",")
        stops = lList[9].split(",")
        symbol = enstSymbol[lList[0]]
        if currentSymbol=="":
            currentSymbol = symbol 
            currentChrom = lList[1]
            currentStrand = lList[2]
            minStart = lList[3]
            maxStop = lList[4]
        #check if the line belongs to the current symbol
        if symbol != currentSymbol :
            
            c.execute("INSERT INTO gene VALUES("+str(geneIdx)+", '""', '"+currentSymbol+"', '"+currentChrom+"', '"+minStart+"', '"+maxStop+"', Null, 'UCSC KnownGene')")
            currentSymbol = symbol 
            currentChrom = lList[1]
            currentStrand = lList[2]
            minStart = lList[3]
            maxStop = lList[4]
            print(geneIdx)
            geneIdx+=1
        #always insert the transcript info
        c.execute("INSERT INTO Transcript VALUES("+str(geneIdx)+", '"+str(tranIdx)+"', '"+lList[0]+"', '"+lList[1]+"', '"+lList[3]+"', '"+lList[4]+"')")
        minStart = minStart if lList[3] > minStart else lList[3]
        maxStop = maxStop if lList[4] < maxStop else lList[4]
            
        """
        c.execute("SELECT id FROM gene WHERE symbol = '"+symbol+"'")
        ret = c.fetchall()
        gid = ret[0][0]
        
        #the number of splices is the length of the starts list -1, but all the lists ends in a comma and artificially inflate the length by 1 so it must be -2 here
        for index in range(len(starts)-2): 
            try:
                c.execute("INSERT INTO splice VALUES("+str(gid)+", '"+str(spliceIdx)+"', 'None', '"+stops[index]+"', 'None', '"+starts[index+1]+"')")
                spliceIdx+=1
            except IntegrityError as e:
                #duplicate splice
                print(e.__str__(), str(gid), stops[index], starts[index+1])
        """
    c.execute("commit")
    

    print('done')
