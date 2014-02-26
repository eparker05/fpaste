#!/home/wwwadmin/fastadb/python27/bin/python
from urllib import quote_plus
import requests

def get_uniprot_protein_info(uniprotId, cutSignalSequence=False):
    uniUrl = "http://www.uniprot.org/uniprot/" + quote_plus(uniprotId) + ".txt"
    result = requests.get(uniUrl)
    if result.status_code != 200:
        return {"success":False, 
                "error":"Http code {} for id {}.".format(result.status_code, uniprotId),
                "header":uniprotId,
                "meta":uniUrl}
    else:
        #transcribe into split lists
        rawResultList = result.text.split("\n")
        uniprotSourceList = []
        for line in rawResultList:
            outputList = line.split()
            uniprotSourceList.append(outputList)
        #find sequence
        foundSQ = False
        for lineList in uniprotSourceList:
            if lineList and lineList[0] == "SQ":
                SqLength = int(lineList[2].strip())
                SqIndex = uniprotSourceList.index(lineList)
                sequence = ''.join([''.join(a) for a in uniprotSourceList[SqIndex+1:-1]])
                sequence = ''.join([a for a in sequence if a.isalpha()])
                finalSequence = sequence
        #No sequence length error state
        if not SqLength:
            return {"success":False, 
                    "error":"Sequence length not found in SQ section",
                    "header":uniprotId,
                    "meta":uniUrl}
        #meta data does not match sequence error state
        if len(sequence) != SqLength:
            return {"success":False, 
                    "error":"Sequence length {} does not match SQ annotation {}".format(len(sequence), SqLength),
                    "header":uniprotId,
                    "meta":uniUrl}
        #find signal sequence when required
        if cutSignalSequence:
            for lineList in uniprotSourceList:
                if lineList and lineList[0] == "FT" and lineList[1] == "SIGNAL":
                    beginSignal = int(lineList[2].strip()) - 1
                    endSignal = int(lineList[3].strip())
                    finalSequence = sequence[endSignal:]        
        return {"sequence":finalSequence, 
                "header":uniprotId, 
                "meta":uniUrl, 
                "success":True, 
                "error":""}
                
_testseq = "P12763"

