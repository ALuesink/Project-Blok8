from flask import Flask, request, render_template
from werkzeug.contrib.cache import SimpleCache
from Bio import Entrez
from Entrez import findArticles, getAbstracts, tabel
from textmining import openbestand, tokenize, overeenkomst, toJson

app = Flask(__name__)
cache = SimpleCache()#dit is een cache om informatie tijdelijk in op te slaan

#Geeft de template index.html, de hoofdpagina, terug.
@app.route('/')
def webintro():
    return render_template('./index.html')

#Deze functie geeft  de template table.html met de tabel met alle resultaten terug
#hierbij wordt het zoekwoord opgehaald, om artikelen te vinden die dit zoekwoord bevatten. 
#Hierbij wordt er gebruik gemaakt van een cache, in deze cache zit een dictionary met belangrijke keywords uit de artikelen.
#Als deze cache leeg, wordt de dictionary_text toegevoegd, die opgehaald is bij getAbstracts(). 
@app.route('/tabel', methods=["GET"])
def tabelweergeven():
    zoekwoord = request.args.get("zoekwoord") 
    Entrez.email = "your.email@example.com"
    webintro() 
    ID = findArticles(zoekwoord)
    keys, abstracts, auteur, datum, titel, dictionary, dictionary_text = getAbstracts(ID) #haalt van elke artikel informatie op en voegt elk soort informatie(bijv. auteurs) in een aparte lijst toe.
    dic = cache.get('dictionary')
    if dic is None:
        dic = dictionary_text
        cache.set('dictionary',dic)
    return render_template('table.html', zoekwoord=zoekwoord, tabel=tabel(ID, auteur, keys, datum, titel))

#Maakt een graph aan met overeenkomstige keywords uit de artikelen. 
#Deze functie maakt gebruik van een cache, waarin een dictionary zit met alle keywords.
#Deze functie geeft de template graph.html terug met de graph.
@app.route('/graph')
def graph():
    dictionary = cache.get('dictionary')
    terrier = openbestand()
    combos = tokenize(dictionary, terrier)
    overeenkomsten = overeenkomst(combos)
    graph_data, json_bestand= toJson(overeenkomsten)
    return render_template('graph.html', data=graph_data, json=json_bestand)


app.secret_key = 'sajfldjslafjlsdajfl;sadjfl;sjaf'

if __name__ == '__main__':
    app.run()
