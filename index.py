from flask import Flask, request, render_template
from werkzeug.contrib.cache import SimpleCache
from Bio import Entrez
from Entrez import findArticles, getAbstracts, tabel
from textmining import openbestand, tokenize, overeenkomst, toJson

app = Flask(__name__)
cache = SimpleCache()

@app.route('/')
def webintro():
    return render_template('./index.html')

@app.route('/tabel', methods=["GET"])
def tabelweergeven():
    zoekwoord = request.args.get("zoekwoord")
    Entrez.email = "your.email@example.com"
    webintro()
    ID = findArticles(zoekwoord)
    keys, abstracts, auteur, datum, titel, dictionary, dictionary_text = getAbstracts(ID)
    dic = cache.get('dictionary')
    if dic is None:
        dic = dictionary_text
        cache.set('dictionary',dic)
    return render_template('table.html', zoekwoord=zoekwoord, tabel=tabel(ID, auteur, keys, datum, titel))

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
