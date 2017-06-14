from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
from nltk.stem import *
from collections import Counter
import string
import re
import json
import os

#functie: openen bestand met extra stopwoorden
#return: lijst met stopwoorden
def openbestand():
    terrierbestand = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'Project_8', 'static', 'stop_terrier.txt'))
    terrier_bestand = open(terrierbestand,"r")

    terrier_stop = []
    line = terrier_bestand.readline()
    while line != "":
        line.rstrip()
        terrier_stop.append(line[:-1])
        line = terrier_bestand.readline()

    terrier_bestand.close()
    return terrier_stop

#input: dictionary van artikelen met titel, abstact en keywoorden, lijst met extra stopwoorden
#functie: per artikel worden de titel, abstract en keywoorden gestemd en bekeken of deze niet in de stopwoorden zit,
#van de 10 meest voorkomende woorden uit de abstract met de titel en keywoorden worden combinaties gemaakt van alle woorden per artikel.
#return: per artikel/PMID combinaties van termen
def tokenize(artikels,terrier_stop):
    leestekens = list(string.punctuation)
    extra_woorden = ["cancer","human","amino","aa","enzym","acid","fatti","vitro","rat","level","death","tumor"]
    stopwoorden = stopwords.words("english") + leestekens + terrier_stop + extra_woorden
    stemmer = SnowballStemmer("english")

    combos = {}
    for key, value in artikels.items():
        titel = value[0]
        abstract = value[1]
        keywords = value[2]
        titelwoorden = []
        abstractwoorden = []
        keywords_woorden = []
        woorden_artikel = []
        combi_artikel = []
        pattern = '[0-9]+(-lipoxygenas)|[0-9]+(-lox)|(lox)|(lipoxygenas)|(cox)|(cyclooxygenas)|(death)|(apoptosi)|(cancer)|(tumor)|(carinoma)'
        synoniemen = {"lox":"lipoxygenas",
                      "lipoxygenas": "lipoxygenas",
                      "-lox": "-lipoxygenas",
                      "cox": "cyclooxygenas",
                      "cyclooxygenas": "cyclooxygenas",
                      "death": "apoptosi",
                      "apoptosi": "apoptosi",
                      "cancer": "cancer",
                      "tumor": "cancer",
                      "carcinoma": "cancer"}

        # print(keywords)
        for term in word_tokenize(titel):
            if term.lower() not in stopwoorden and not term.isdigit():
                titelwoorden.append(stemmer.stem(term.lower()))
        for term in word_tokenize(abstract):
            if term.lower() not in stopwoorden and not term.isdigit():
                abstractwoorden.append(stemmer.stem(term.lower()))
        for term in keywords:
            if term.lower() not in stopwoorden and not term.isdigit():
                keywords_woorden.append(stemmer.stem(term.lower()))

        count_artikel = Counter()
        count_artikel.update(abstractwoorden)
        common10_artikel = count_artikel.most_common(10)

        for i in range(len(common10_artikel)):
            woorden_artikel.append(common10_artikel[i][0])

        woorden_artikel.extend(titelwoorden)
        woorden_artikel.extend(keywords_woorden)


        for i in range(len(woorden_artikel)-1):
            for j in range(i+1, len(woorden_artikel)):
                w1,w2 = sorted([str(woorden_artikel[i]),str(woorden_artikel[j])])
                matchw1 = re.match(pattern,w1,flags=0)
                matchw2 = re.match(pattern,w2,flags=0)
                if matchw1 is not None:
                    w1_match = matchw1.group()
                    positie = w1_match.find("-lox")
                    if positie != -1:
                        aantal = w1_match[:(len(w1_match)-len("-lox"))]
                        w1 = aantal + "-lipoxygenas"
                    else:
                        if bool(re.search(r'[0-9]+',w1_match)):
                            w1 = w1_match
                        else:
                            w1 = synoniemen.get(w1_match)
                if matchw2 is not None:
                    w2_match = matchw2.group()
                    positie = w2_match.find("-lox")
                    if positie != -1:
                        aantal = w2_match[:(len(w2_match)-len("-lox"))]
                        w2 = aantal + "-lipoxygenas"
                    else:
                        if bool(re.search(r'[0-9]+', w2_match)):
                            w2 = w2_match
                        else:
                            w2 = synoniemen.get(w2_match)
                if w1 != w2:
                    combi_artikel.append((w1,w2))

        combos[key] = combi_artikel
    return combos

#input: dictionary van combinaties per artikel/PMID
#functie: vergelijken van de combinaties tussen de artikelen. 
#Alleen de combinaties die in 10 of meer artikelen staan worden opgeslagen in een dictionary
#met als key de combinatie en als values de PMIDs
#return: dictionary van combinaties met PMIDs
def overeenkomst(combos):
    keyList = sorted(combos.keys())
    overeenkomst_dic = {}
    for i in range(len(keyList) - 1):
        for j in range(i + 1, len(keyList)):
            lijst_1 = combos.get(keyList[i])
            lijst_2 = combos.get(keyList[j])

            overeenkomst = list(set(lijst_1).intersection(lijst_2))

            al_in_dic = sorted(overeenkomst_dic.keys())
            for item in range(len(overeenkomst)):
                combi = overeenkomst[item]
                if combi in al_in_dic:
                    lijstIDs = set(overeenkomst_dic.get(combi))
                    lijstIDs.add(keyList[i])
                    lijstIDs.add(keyList[j])
                    overeenkomst_dic[combi] = list(lijstIDs)
                elif combi not in al_in_dic:
                    overeenkomst_dic[combi] = [keyList[i],keyList[j]]

    dic_keys = sorted(overeenkomst_dic.keys())
    eind_dic = {}

    for i in range(len(dic_keys)):
        combi = dic_keys[i]
        values = overeenkomst_dic.get(combi)
        if len(values) >= 10:
            eind_dic[combi] = values

    return eind_dic

#input: dictionary van combinaties met PMIDs
#functie: het omzetten van de input naar JSON format, ook wordt een JSON bestand opgeslagen
#return:2 lijsten met JSON data, nodes en links voor de visualisatie van de graph
def toJson(dic_overeenkomst):
    project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'Project_8', 'static'))
    json_bestand = {"nodes": [], "links": []}

    combis = sorted(dic_overeenkomst.keys())

    set_nodes = set()
    for i in range(len(combis)):
        node_1 = combis[i][0]
        node_2 = combis[i][1]
        set_nodes.add(node_1)
        set_nodes.add(node_2)

    list_nodes = sorted(set_nodes)

    label_id = {}
    for i in range(len(list_nodes)):
        json_bestand["nodes"].append({
            "id": str(list_nodes[i]),
            "loaded": True,
            "name" : str(list_nodes[i]),
            "style": {"label":str(list_nodes[i])}
        })
        label_id[str(list_nodes[i])] = str(list_nodes[i])

    for j in range(len(combis)):
        json_bestand["links"].append({
            "id": str(combis[j]),
            "from": label_id.get(combis[j][0]),
            "to": label_id.get(combis[j][1]),
            "style": {"label": len(dic_overeenkomst.get(combis[j]))}
        })

    with open(os.path.join(project_path, "data.json"),"w") as out_file:
        json.dump(json_bestand,out_file)

    nodes = json_bestand["nodes"]
    links = json_bestand["links"]
    graph_data = [nodes, links]

    return graph_data, json_bestand
