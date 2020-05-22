import os
from datetime import date,datetime,timedelta
from utils import Utils
from elsapy.elsclient import ElsClient
from elsapy.elsprofile import ElsAuthor, ElsAffil
from elsapy.elsdoc import FullDoc, AbsDoc
from elsapy.elssearch import ElsSearch
from Bio import Entrez
import json
from visualization import Visualization

class AcceptanceAnalysis:

    def __init__(self):
        ## Load configuration
        con_file = open("config.json")
        config = json.load(con_file)
        con_file.close()

        ## Initialize client
        self.client = ElsClient(config['apikey'])
        self.client.inst_token = config['insttoken']

    def get_dates(self, file_name):
        contents = utils.read_from_csv(file_name).split('\n')
        contents = contents[0:-1]
        acceptance_rate=dict()
        acceptance_rate['week']=0
        acceptance_rate['fortnight']=0
        acceptance_rate['triweek']=0
        acceptance_rate['month_plus']=0
        acceptance_by_journal=dict()

        unknowns = []
        num_handled_articles = 0
        count_invalid=0
        for idx, line in enumerate(contents):
            if len(line) == 1:
                continue
            title, authors, doi, journal_name, date_online, url = line.split(',')
            res = self.get_data_from_doi(doi)
            # print(res)
            if not isinstance(res,str):
                if not isinstance(res['PubmedArticle'], list):
                    continue
                try:
                    num_handled_articles += 1
                    if not journal_name in acceptance_by_journal.keys():
                        acceptance_by_journal[journal_name]=dict()
                        acceptance_by_journal[journal_name]['total'] =0
                        acceptance_by_journal[journal_name]['num_articles']=0


                    dates = res['PubmedArticle'][0]['PubmedData']['History']
                    for article_date in dates:
                        if article_date.attributes['PubStatus']=='received':
                            date_received=date(int(article_date['Year']), int(article_date['Month']), int(article_date['Day']))
                        if article_date.attributes['PubStatus']=='revised':
                            date_revised = date(int(article_date['Year']), int(article_date['Month']), int(article_date['Day']))
                        if article_date.attributes['PubStatus']=='accepted':
                            date_accepted = date(int(article_date['Year']), int(article_date['Month']), int(article_date['Day']))

                    acceptance_time=date_accepted-date_received+timedelta(days=1)
                    print(acceptance_time)
                    if acceptance_time.total_seconds()<=0:
                        count_invalid+=1
                        continue
                    if acceptance_time.days<=7:
                        acceptance_rate['week']+=1
                    else:
                        if acceptance_time.days<=14:
                            acceptance_rate['fortnight']+=1
                        else:
                            if acceptance_time.days <= 21:
                                acceptance_rate['triweek'] += 1
                            else:
                                acceptance_rate['month_plus'] += 1
                    acceptance_by_journal[journal_name]['total']+=acceptance_time.days
                    acceptance_by_journal[journal_name]['num_articles']+=1

                except :
                    print('caught exception for {}'.format(res))
        print('num handled articles {}, out of {}. count of invalid acceptance date data - {}'.format(num_handled_articles,idx+1,count_invalid))
        for k,v in acceptance_by_journal.items():
            if v['num_articles']==0:
                v['num_articles']=1
            avg=v['total']/v['num_articles']
            print('for journal {}, num articles {}, avg acceptance time in days {}'.format(k,v['num_articles'],avg))
        return acceptance_rate, idx+1,num_handled_articles,count_invalid



    def get_data_from_doi(self, doi):
        doi_doc = FullDoc(doi=doi)
        if doi_doc.read(self.client):
            # print("doi_doc.title: ", doi_doc.title)
            doi_doc.write()
        else:
            print("Read document failed.")
            return doi
        id = None
        if not 'pubmed-id' in doi_doc._data.keys():
            print("no pubmed-id, trying with title")
            # try with title
            Entrez.email = 'shirAviv@protonmail.com'
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
            Entrez.email = 'shirAviv@protonmail.com'
            handle = Entrez.efetch(db='pubmed',
                                   retmode='xml',
                                   id=id)
            results = Entrez.read(handle)
            # print(results)
            return results

        else:
            print("no pubmed id")
            return doi


if __name__ == '__main__':
    print(datetime.now())
    utils = Utils(path='D:\\shir\\study\\covid_19\\scopus')
    aa = AcceptanceAnalysis()
    vis=Visualization()
    data_dict=dict()
    months=['Apr']
    for month in months:
        data_dict[month]=dict()
        acc_speed, num_total_articles, num_handled_articles, num_invalid_speed_articles=aa.get_dates(month+'_articles.csv')
        print(acc_speed)
        data_dict[month]['Total']=num_handled_articles-num_invalid_speed_articles
        data_dict[month]['acc_speed']=acc_speed
    # vis.plt_acceptance_speed(data_dict,"acceptance speed")
    print(datetime.now())
