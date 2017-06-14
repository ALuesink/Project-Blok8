from Bio import Entrez, Medline

#De zoekopdracht in pubmed. Er wordt een connectie gemaakt met behulp van Entrez.
#Er is ingesteld dat alle artikelen vanaf het jaar 2000 worden gezocht. 
def findArticles(zoekwoord):
    handle = Entrez.esearch(db="pubmed",
                            term=zoekwoord,
                            mindate='2000',
                            usehistory='y',
                            retmax='100000')
    record = Entrez.read(handle)
    handle.close()
    ID = record["IdList"]

    return ID

#Met behulp van de efetch functie van Entrez zijn de abstracts, keywords, auteurs, data en titels van elk artikel opgehaald. 
def getAbstracts(ID):
    dictionary = {}
    dictionary_textmining = {}
    abstracts = []
    keys = []
    auteur = []
    datum = []
    titel = []

    handle = Entrez.efetch(db="pubmed",
                           id=ID,
                           rettype='Medline',
                           retmode='text')
    records = Medline.parse(handle)
    for record in records:
        PMID = record.get('PMID')
        auteurs = record.get('AU')
        if record.get('AB') is not None:
            abstract = record.get('AB')
        else:
            abstract = "-"
        date = record.get('DP')
        title = record.get('TI')
        if record.get('OT') is None:
            keywords = "-"
        else:
            keywords = record.get('OT')

        auteur.append(auteurs)
        abstracts.append(abstract)
        datum.append(date)
        titel.append(title)
        keys.append(keywords)

        dictionary[PMID] = [title,abstract,keywords,auteurs,date]
        dictionary_textmining[PMID] = [title,abstract,keywords]

    return keys, abstracts, auteur, datum, titel, dictionary, dictionary_textmining

#Generen van de tabel met de ID's, auteurs, keywords, data en titels. 
def tabel(ID, auteur, keys, datum, titel):
    tablet = """"""
    aantal = 0
    for item in ID:
        url = ("http://www.ncbi.nlm.nih.gov/pubmed/%s" % str(item))
        regel = "<tbody><tr><td><font color=""black""><a href=" + url + " target=""_blank""></font>" + str(
            item) + "</a></td><td>" + str(titel[aantal]) + "</td><td>" + ", ".join(auteur[aantal]) + "</td><td>" + str(
            datum[aantal]) + "</td><td><div class=""comment more"">" + ", ".join(keys[aantal]) + "</div></td></tbody>"
        tablet += regel
        aantal += 1
    return tablet
