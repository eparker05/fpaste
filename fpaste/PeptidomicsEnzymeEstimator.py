"""
This work is licensed under a Creative Commons Attribution 3.0 Unported License.
http://creativecommons.org/licenses/by/3.0/
"""

from Bio import SeqIO
import re


_el = { "Arg-C proteinase":                r'([A-Z_]R[A-Z][A-Z_])',
        "Asp-N endopeptidase":             r'([A-Z_][A-Z]D[A-Z_])',
        "BNPS-Skatole":                    r'([A-Z_]W[A-Z][A-Z_])',
        "Chymotrypsin low-specificity":    r'([A-Z_][LYF][A-OQ-Z][A-Z_]|[A-Z_]W[A-LNOQ-Z][A-Z_]|[A-Z_]M[A-OQ-W][A-Z_]|[A-Z_]H[B-OQ-W][A-Z_])',
        "Chymotrypsin high-specificity":   r'([A-Z_][YF][A-OQ-Z][A-Z_]|[A-Z_]W[A-LNOQ-Z][A-Z_])',
        "Pepsin":                          r'([A-Z_][FWYL][A-Z][A-Z_]|[A-Z_][A-Z][FWYL][A-Z_])',
        "Plasmin":                         r'([A-Z_][KR][A-Z][A-Z_])',
        "Cathepsin D":                     r'([A-Z_][AVLIPMFW][AVLIPMF][A-Z_])',
        "Thrombin":                        r'([A-Z_]RG[A-Z_]|GR[A-Z][A-Z_])',
        "Elastase":                        r'([A-Z_][AVLIGR][GPALF][A-Z_])',
        "Trypsin":                         r'([A-Z_][KR][A-OQ-Z][A-Z_])'}
        
EnzymeList = {enz: re.compile(_el[enz]) for enz in _el}
enzymeListNC = {enz: {"n":re.compile(r'\A'+_el[enz]), "c":re.compile(_el[enz]+r'\Z')} for enz in _el}


class Peptide(object):
    """
    This class represents the peptide, contains macroscopic information and search functions
    
    Attributes:
    .sequence = a string
    .accession = a string representing a protein accession
    .contextSequence = a string where the leading and trailing two AAs are represented
    .intensity = float of intensity
    .rt = float of retention time (defaults to 0)
    .nMatches = []; names of matched regex patterns on n-side cleavages
    .cMatches = []; names of matched regex patterns on c-side cleavages
    .start = int of start index in original protein
    .end = int of end index in original 

    Methods:     
    .__init__(sequence="", intensity=1, rt=0, sampleId=None, accession=None, proteinDict=[], enzymeRegexDict=[])
        Only use enzymeRegexDict and proteinList togeather 
        
    .assign_protein(proteinDictList):
        proteinDictList = {"accession":"SEQUENCE", ...}
        This method will assign .accession .contextSequence .start and .end
    
    .assign_cleavages(regexDict)
        regexDict = {"enzymeName":"XYZW", ...} where XYZW is a four character
                                pattern with cleavage between Y and Z
        This method will assign .nMatches and .cMatches from the dictionary
    """
    
    def __init__(self, sequence="", intensity=1, rt=0, sampleId=None, accession=None, proteinDict={}, enzymeRegexDictNC={}):
        self.sequence = sequence
        self.intensity = intensity
        self.rt = rt
        self.sampleId = sampleId
        self.accession = accession
        
        self.contextSequence = None
        self.start = None
        self.end = None
        self.assign_protein(proteinDict)
        
        self.nMatches = []
        self.cMatches = []
        if self.contextSequence is not None:
            self.assign_cleavages(enzymeRegexDictNC)
            
    def assign_protein(self, proteinDict):           
        peptideRe = re.compile(self.sequence) 
        #shorten the search space if there is a valid accession id
        if self.accession is not None and self.accession in proteinDict:
            newSingleProteinDict = {self.accession: proteinDict[self.accession]}
        else:
            if self.accession is None:
                raise LookupError("Accession must be assigned for peptide {}".format())
            else:
                raise LookupError("Accession '{}' was not included in the library".format(self.accession))
        #search the protein for the sequence, generate the context sequence
        for proteinKey in newSingleProteinDict:
            sequence = proteinDict[proteinKey]
            location = peptideRe.search(sequence)
            if location is not None:
                self.start = location.start()
                self.end = location.end()
                #figure out the padding situation
                tempStart = max(0, location.start()-2)
                tempEnd = min(len(sequence), location.end()+2)
                startPadding = tempStart - (location.start() - 2)
                endPadding = tempEnd - (location.end() + 2)
                self.contextSequence = startPadding*"_" + sequence[tempStart:tempEnd] + endPadding*"_"
                self.accession = proteinKey
                return True
        #on failure to find sequence, return false
        return False
    
    def assign_cleavages(self, enzymeRegexDictNC, greedy=True):
        #fail immediately if no context sequence exists
        if self.contextSequence is None:
            return False
        #search through regex strings for nSide Matches
        contextSequence = self.contextSequence
        for enzyme in enzymeRegexDictNC:
            cSideRe = enzymeRegexDictNC[enzyme]["c"]
            nSideRe = enzymeRegexDictNC[enzyme]["n"]
            if cSideRe.search(contextSequence) is not None:
                self.cMatches.append(enzyme)
            if nSideRe.search(contextSequence) is not None:
                self.nMatches.append(enzyme)
        
    def __repr__(self):
        return "<Peptide: accession={}, sampleId={} range=({},{}) >".format(self.accession, self.sampleId, self.start, self.end)

        
 
 
def import_peptides_and_preprocess(peptideCsvFileObject, fastaFileObject, validEnzymeList=[]):
    """
    This function takes a format compliant CSV file object, a fasta file object
    and an enzyme regex dictionary containing n side and c side varriations
    """
    if len(validEnzymeList) > 0:
        enzymeRegexDictNC = {enzyme:enzymeListNC[enzyme] for enzyme in validEnzymeList if enzyme in enzymeListNC}
    else:
        enzymeRegexDictNC = enzymeListNC
    proteinDict = {p.id:str(p.seq) for p in SeqIO.parse(fastaFileObject, "fasta")}
    firstLinePending = True
    peptideList = []
    
    headerDict = {"sequence":False, "intensity":False, "protein_id":False, "sample_id":False, "rt":False}
    fileAnalysisList = {}
    for line in peptideCsvFileObject:
        #parse header indicies
        if firstLinePending and (len(line) > 35 and "#" not in line[0:3]):
            firstLinePending = False
            splitLine = [li.lower().strip() for li in line.split(',')]
            for i in range(len(splitLine)):
                if splitLine[i] in headerDict:
                    headerDict[splitLine[i]] = i
            if (not headerDict["intensity"] and not headerDict["sequence"]
                    and not headerDict["protein_info"] and not headerDict["sample_id"]):
                raise ValueError("Incorrect file format; missing one or more columns")
            continue
        #try extracting the data, failures are silent and will result in the line/peptide being ignored        
        try:
            splitLine = [li.strip() for li in line.split(',')]
            rt = None
            if headerDict["rt"]:
                rt = float(splitLine[headerDict['rt']])
            sequence = splitLine[headerDict['sequence']]
            intensity = float(splitLine[headerDict['intensity']])
            sampleId = splitLine[headerDict['sample_id']]
            proteinId = splitLine[headerDict['protein_id']]
        except:
            print "what?"
            continue
        peptide = Peptide(sequence = sequence,
                          intensity = intensity,
                          rt = rt,
                          sampleId = sampleId,
                          accession = proteinId,
                          proteinDict = proteinDict,
                          enzymeRegexDictNC = enzymeRegexDictNC)
        peptideList.append(peptide)
    
    if len(peptideList) == 0:
        raise ValueError("No peptides found, check file format")
    return peptideList

    
def extract_data_from_processed_peptides(peptideList, validEnzymeList, method="sampleId", result="dictionary"): # extractOrderSet,
    """
    This function takes a format compliant CSV file object, a fasta file object
    and an enzyme regex dictionary containing n side and c side varriations
    """
    if len(validEnzymeList) > 0:
        enzymeDict = {enzyme:enzymeListNC[enzyme] for enzyme in validEnzymeList if enzyme in enzymeListNC}
    else:
        enzymeDict = enzymeListNC
        
    outDict = {}
    aaList = [ "A","G","P","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T" ]                          
    
    for peptide in peptideList:
        attribute = getattr(peptide,method)
        
        if attribute not in outDict:
            outDict[attribute] = {"cSideOrphans":{ aa:0 for aa in aaList },
                                  "nSideOrphans":{ aa:0 for aa in aaList },
                                  "enzymeResponseDict":{ enzyme: 0 for enzyme in enzymeDict } }
        
        for enzyme in peptide.nMatches:
            outDict[attribute]["enzymeResponseDict"][enzyme] += peptide.intensity
        for enzyme in peptide.cMatches:
            outDict[attribute]["enzymeResponseDict"][enzyme] += peptide.intensity
        
        if len(peptide.nMatches) < 1:
            if peptide.contextSequence[1] != '_':
                outDict[attribute]["nSideOrphans"][peptide.contextSequence[2]] += 0.5*peptide.intensity
            if peptide.contextSequence[-2] != '_':
                outDict[attribute]["nSideOrphans"][peptide.contextSequence[-2]] += 0.5*peptide.intensity
        
        if len(peptide.cMatches) < 1:
            if peptide.contextSequence[-2] != '_':
                outDict[attribute]["cSideOrphans"][peptide.contextSequence[-3]] += 0.5*peptide.intensity
            if peptide.contextSequence[1] != '_':
                outDict[attribute]["cSideOrphans"][peptide.contextSequence[1]] += 0.5*peptide.intensity
            
    if result == "dictionary":
        return outDict
    
    validAttributes = [attr for attr in outDict]
    validEnzymes = [enz for enz in enzymeDict]
    
    firstColumn = [""] + validEnzymes + [aa+"-c-side" for aa in aaList] + [aa+"-n-side" for aa in aaList]
    columnList = [firstColumn]
    for attr in validAttributes:
        enzResponses = outDict[attr]["enzymeResponseDict"]
        cSide = outDict[attr]["cSideOrphans"]
        nSide = outDict[attr]["nSideOrphans"]
        
        nextColumn = [attr] + [enzResponses[enz] for enz in enzResponses]
        nextColumn += [cSide[enz] for enz in cSide]
        nextColumn += [nSide[enz] for enz in nSide]
        columnList.append(nextColumn)
    
    newList = zip(*columnList)
    newList = [list(tup) for tup in newList]
    if result == "list":
        return newList






"""if __name__ == "__main__":
    import time
    
    t1 = time.time()
    #present to user:
    enzymeList = [enzyme for enzyme in enzymeListNC]
    
    #User input:
    inputEnzymeList = ["Elastase", "Plasmin"] #any valid enzymes from enzymeList
    testFile = "test_15mothers_degradome - Copy.csv"
    #testFile = "test_15mothers_degradome.csv"
    testFile = open(testFile, 'r')
    
    
    #prepare output after successfull request
    fastaLibraryFile = open("degradomeLibrary.fasta", "r")
    #peptideList = import_peptides_and_preprocess(testFile, fastaLibraryFile, inputEnzymeList)
    try:
        peptideList = import_peptides_and_preprocess(testFile, fastaLibraryFile, inputEnzymeList)
    except Exception, e:
        print e
    try:
        results = extract_data_from_processed_peptides(peptideList, enzymeListNC, result="list")
    except Exception, e:
        print e
        
    print "took {} seconds".format(time.time() - t1)"""


if __name__ == "__main__":
    #from EnzymeCutQuantifier import extract_data_from_processed_peptides, import_peptides_and_preprocess
    
    fastaFile = open("degradomeLibrary.fasta", "r")
    csvTestFile = open("test_15mothers_degradome.csv", 'r')
    inputEnzymeList = ["Trypsin", "Elastase"]
    peptideList = import_peptides_and_preprocess(csvTestFile, fastaFile, inputEnzymeList)
    resultCSV = extract_data_from_processed_peptides(peptideList, inputEnzymeList, result = "list")
    
    import csv
    import csv
    with open('!output.csv', 'wb') as csvfile:
        resultWriter = csv.writer(csvfile)
        resultWriter.writerows(resultCSV)