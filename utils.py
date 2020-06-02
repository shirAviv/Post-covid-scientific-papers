import os
import pickle
from elsapy.elsclient import ElsClient
from elsapy.elsprofile import ElsAuthor, ElsAffil
from elsapy.elsdoc import FullDoc, AbsDoc
from elsapy.elssearch import ElsSearch
from Bio import Entrez
import json
from pybliometrics.scopus import ScopusSearch



class Utils():
    def __init__(self,path):
        self.path=path
        con_file = open("config.json")
        config = json.load(con_file)
        con_file.close()

        ## Initialize client
        self.client = ElsClient(config['apikey'])
        self.client.inst_token = config['insttoken']

    def save_obj(self,obj, name):
        obj_path = os.path.join(self.path, name)
        # str=pickle.dumps(obj)
        with open(obj_path + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.DEFAULT_PROTOCOL)

    def load_obj(self,name):
        obj_path = os.path.join(self.path, name)
        with open(obj_path + '.pkl', 'rb') as f:
            return pickle.load(f)

    def write_to_csv(self, df, name):
        df.to_csv(name, index=False)

    def read_from_csv(self,name):
        csv_path=os.path.join(self.path, name)
        with open(csv_path, encoding="utf8") as csv_file:
            contents = csv_file.read()
        return contents

    def get_data_from_doi(self, doi,title):
        id=None
        affil=None
        try:
            doi_doc=ScopusSearch(doi, subscriber=False)
            if 'pubmed-id' in doi_doc._json[0].keys():
                id=doi_doc._json[0]["pubmed-id"]
            if 'affiliation' in doi_doc._json[0].keys():
               affil= doi_doc._json[0]['affiliation']
            if id==None:
                doi_doc = FullDoc(doi=doi)
                if doi_doc.read(self.client):
                    # print("doi_doc.title: ", doi_doc.title)
                    doi_doc.write()
                else:
                    print("Read document failed. no id for doi {}. trying with title".format(doi))
                    doi_doc=None
                    # return doi, affil
                id = None
                if doi_doc==None or (not 'pubmed-id' in doi_doc._data.keys()):
                    print("trying with title")
                    # try with title
                    Entrez.email = 'shirAviv@protonmail.com'
                    if doi_doc==None:
                        query = title
                    else:
                        query = doi_doc.title
                    handle = Entrez.esearch(db='pubmed',
                                            retmode='xml',
                                            term=query)
                    results = Entrez.read(handle)
                    if int(results['Count']) > 0:
                        id = results['IdList']
                else:
                    id = doi_doc._data['pubmed-id']
            if id != None:
                return self.fetch_data_from_pubmed(id), affil

            else:
                print("no pubmed id found for doi {}".format(doi))
                return doi, affil
        except:
            print("caught exception for doi {}".format(doi))
            return doi, affil

    def fetch_data_from_pubmed(self, id):
        Entrez.email = 'shirAviv@protonmail.com'
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=id)
        results = Entrez.read(handle)
        # print(results)
        return results
