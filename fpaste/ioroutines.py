#!/home/wwwadmin/fastadb/python27/bin/python
from urllib import quote_plus
import requests
from Bio import SeqIO
from io import StringIO

def get_uniprot_protein_info(uniprotId, cutSignalSequence=False):
    uniUrl = "http://www.uniprot.org/uniprot/" + quote_plus(uniprotId) + ".txt"
    page = requests.get(uniUrl)
    outDict = { "success":False, 
                "error":"",
                "header":uniprotId,
                "meta":uniUrl }
    #check for bad return
    if page.status_code != 200:
        outDict["error"] = "Http code {} for id {}.".format(page.status_code, uniprotId)
        return outDict
    #create a file object for the parser
    uniFileFake = StringIO(page.text)
    try:
        seqRec = SeqIO.read(uniFileFake, 'swiss')
    except ValueError:
        outDict["error"] = "Invalid format was returned for id {}.".format(uniprotId)
        return outDict
    #parse output for SIGNAL sequence, be very paranoid about failure, return 0 by default or in failure 
    realStart = 0
    if cutSignalSequence:
        for feat in seqRec.features:
            if feat.type == "SIGNAL":
                try:
                    realStart = feat.location.end
                except:
                    pass
    sequence = str(seqRec.seq)[realStart:]
    outDict["seq"] = sequence
    outDict["success"] = True
    return outDict