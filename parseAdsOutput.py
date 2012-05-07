from string import *
#------------------------------------------------------------------
# generates a csv file containing info about the publication
# and the names and the authors who are from leiden
# there is a problem getting the individual authors!! for example when 
# an autor is listen only as last name i.e not in the form "last,first" 
# but as "author," the order is getting messed up
#------------------------------------------------------------------
origNameListPath = '/home/mher/adsStuff/nameList2011Mher.txt'

#fName        = '/home/mher/adsStuff/sample.ads'
#outCsvFname  = '/home/mher/adsStuff/sample.csv'

fName        = '/home/mher/adsStuff/non-refereed.ads'
outCsvFname  = '/home/mher/adsStuff/non-refereed.csv'

#fName        = '/home/mher/adsStuff/refereed.ads'
#outCsvFname  = '/home/mher/adsStuff/refereed.csv'


pubs         = ()  # holds all the publication items
inputAuthors = ()  # holds all the names which were used as input for the ads search 
#----------------------------------------------------------------
class author():
    def __init__(self):
        self.rank  = None
        self.first = None
        self.Last  = None
#_---------------------------------------------------------------
class publication():
    
    def __init__(self):
        self.title    = None
        self.authors  = None
        self.date     = None
        self.journal  = None
        self.volume   = None
        self.page     = None
        
        self.nAuth    = None
        self.authsObj = None
        
    def parseFromString(self, strng):

        def parseAuthors( authStrng ):
            
            authsObj = ()
            authStrng = authStrng.replace(", and ", ",") #if this is not found, strng doest change
            authStrng = authStrng.replace(" and ", ",")  #if above is found, this will not be found anyway
            strngSplt = authStrng.split(',')
            
            nAuth = 0
            nParts = len(strngSplt)     
            for i, cmpnt in enumerate(strngSplt):
                
                if i % 2 == 0:
                    auth = author()
                    
                     
                    #print i,
                    auth.rank = (i/2)  + 1
                    
                    #print strngSplt[i].strip(),
                    auth.last = strngSplt[i].strip()
                    
                    if (i+1) < nParts:
                        #print strngSplt[i+1].strip()
                        auth.first = strngSplt[i+1].strip()
                    else:
                        auth.first = " "
                        #print
                    nAuth += 1
                    
                    authsObj += (auth,)
                else:
                    continue
            
            return (nAuth, authsObj) 
            
                
        strng = strng.replace("\n", " ")
        itemStrCmpnt = strng.split("XXXXX")
        if len(itemStrCmpnt) == 6 :
            # getting the components of the item        
            self.title   = itemStrCmpnt[0].strip()
            self.authors = itemStrCmpnt[1].strip()
            self.date    = itemStrCmpnt[2].strip()
            self.journal = itemStrCmpnt[3].strip()
            self.volume  = itemStrCmpnt[4].strip()
            self.page    = itemStrCmpnt[5].strip()
            # parsing the author names
            self.nAuth, self.authsObj, = parseAuthors(self.authors)
        #print self.title
    
    def show(self):
        print 'title    = %s\nAuthors = %s\nDate     = %s\nJournal  = %s\nVolume   = %s\npage     = %s\n' % (self.title, self.authors, self.date, self.journal, self.volume, self.page)
#-------------------------------------------------------------------------------
def writeCSV(pubs, outFileName):
    
    fObj = open(outFileName, 'w')
    
    for item in pubs:

        #                0   1    2     3   4    5    6    7    8     9  
        strngFormat = '"%s","%s","%s","%s","%s","%s","%d","%s","%s","%d"\n'
        strngContnt = ()                                      
        strngContnt += ( (item.title).replace('"', ' '),)     # 0
        strngContnt += (item.authors,)                        # 1
        strngContnt += (item.date,)                           # 2
        strngContnt += (item.journal,)                        # 3
        strngContnt += (item.volume,)                         # 4
        strngContnt += (item.page,)                           # 5
                  
        strngContnt += (item.nAuth,)                          # 6
        
        
        """
        # checking if leiden authors are among the authors
        # this is the faster but less accurate version of the code below
        leidenAuthors  = ()
        firstAuthPaper = 0
        for i, auth1 in enumerate(item.authsObj):
            for auth2 in inputAuthors:
                if auth1.last == auth2.last:
                    leidenAuthors += (auth2,)
                    if i == 0:
                        firstAuthPaper = 1
        """
  
        nLeidenAuth = 0
        leidenAuthStrng = ""        
        for leidenAuth in inputAuthors:
            if leidenAuth.last in item.authors:
                leidenAuthStrng += leidenAuth.last + ',' + leidenAuth.first + ','
                nLeidenAuth += 1
        
        strngContnt += (item.authsObj[0].last,) # 7     
        strngContnt += (leidenAuthStrng,)       # 8
        strngContnt += (nLeidenAuth,)           # 9
        

        strng =  strngFormat % strngContnt
        
        fObj.write(strng) 
    
    fObj.close()
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# ************ reading and parsing the ads output ********************
# reading the whole file into a single big string
# while skipping the header which is 5 lines
fObj = open(fName, 'r')
fileStr = ''
lineNum = 0
for line in fObj:
    if lineNum >= 5: 
        fileStr += line
    lineNum += 1
fObj.close()

# parsing each item into a tuple
#items print fileStr
itemsStr = fileStr.split('NNNNN')

# converting the string items to parsed publication objects
for itemStr in itemsStr:
    item = publication()
    item.parseFromString(itemStr)
    if item.title != None:
        pubs = pubs + (item,) 

# ************ reading and parsing the input name list ********************
fObj = open(origNameListPath, 'r')

for line in fObj:
    auth = author()
    
    lineSplt = line.split(',')
    auth.first = (lineSplt[1].replace('\n', ' ')).strip() 
    auth.last  = lineSplt[0]
    
    inputAuthors += (auth,) 

# ************ reading and parsing the input name list ********************
writeCSV(pubs, outCsvFname)

print len(pubs)
for auth in inputAuthors:
    print auth.last, '---', auth.first
    
for item in pubs:
    item.show() 
print 'done'