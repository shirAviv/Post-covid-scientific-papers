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
                        acceptance_by_journal[journal_name]['total_days'] =0
                        acceptance_by_journal[journal_name]['num_articles']=0
                        acceptance_by_journal[journal_name]['ids']=set()

                    id = res['PubmedArticle'][0]['MedlineCitation']['PMID']
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
                    acceptance_by_journal[journal_name]['total_days']+=acceptance_time.days
                    acceptance_by_journal[journal_name]['num_articles']+=1
                    acceptance_by_journal[journal_name]['ids'].add(id)

                except :
                    print('caught exception for {}'.format(res))
        print('num handled articles {}, out of {}. count of invalid acceptance date data - {}'.format(num_handled_articles,idx+1,count_invalid))
        for k,v in acceptance_by_journal.items():
            if v['num_articles']==0:
                v['num_articles']=1
            avg=v['total_days']/v['num_articles']
            print('for journal {}, num articles {}, avg acceptance time in days {}'.format(k,v['num_articles'],avg))
        return acceptance_rate, acceptance_by_journal, idx+1,num_handled_articles,count_invalid

    def fetch_data_from_pubmed(self,id):
        Entrez.email = 'shirAviv@protonmail.com'
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=id)
        results = Entrez.read(handle)
        # print(results)
        return results

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
           return self.fetch_data_from_pubmed(id)

        else:
            print("no pubmed id")
            return doi

    def get_journal_data(self, journal_name, min_date, max_date):
        query = journal_name + '[journal]'
        # query='Travel+Medicine+and+Infectious+Disease[journal]'
        Entrez.email = 'shirAviv@protonmail.com'
        handle = Entrez.esearch(db='pubmed',
                                sort='date',
                                retmax='200',
                                retmode='xml',
                                mindate=min_date,
                                maxdate=max_date,
                                term=query)
        results = Entrez.read(handle)
        return results


    def iterate_journals(self, covid_data_dict):
        acceptance_by_journal = dict()
        empty_results=set(())
        journals=utils.read_from_csv('Journals_by_month_keyword_short_2.csv').split('\n')
        journals=journals[0:-1]
        for idx,journal in enumerate(journals):
            if idx==10:
                break
            journal_name,value,articles_in_press=journal.split(',')
            if journal_name.lower().endswith('jan'):
                year = 2020
                month = 1
                month_str='Jan'
                continue
            if journal_name.lower().endswith('feb'):
                year = 2020
                month = 2
                month_str='Feb'
                continue
            if journal_name.lower().endswith('mar'):
                year = 2020
                month = 3
                month_str='Mar'
                continue
            if journal_name.lower().endswith('apr'):
                year = 2020
                month = 4
                month_str='Apr'
                continue
            if not journal_name in acceptance_by_journal.keys():
                acceptance_by_journal[journal_name]=dict()
            str_year_month = str(year) + '/' + str(month)
            acceptance_by_journal[journal_name][str_year_month]=dict()
            acceptance_by_journal[journal_name][str_year_month]['total_days'] = 0
            acceptance_by_journal[journal_name][str_year_month]['num_articles'] = 0

            if journal_name.lower().endswith('the lancet'):
                empty_results.add(journal_name)
                continue
            res_journals = self.get_journal_data(journal_name, str_year_month, str_year_month)
            # num_publication[journal_name][str_year_month]['total'] = str(res_journals['Count'])
            if int(res_journals['Count'])==0:
                tmp_journal_name=journal_name.replace("&","")
                res_journals = self.get_journal_data(tmp_journal_name, str_year_month, str_year_month)
                # num_publication[journal_name][str_year_month]['total'] = str(res_journals['Count'])
                if int(res_journals['Count']) == 0:
                    empty_results.add(journal_name+" "+str(month))

            print('working on {}'.format(journal_name))
            covid_acceptance_data=data_dict[month_str]['journals'][journal_name]
            acceptance_by_journal=self.get_paper_dates(acceptance_by_journal, journal_name, res_journals, str_year_month, covid_acceptance_data)


            year_older = 2019
            str_year_month = str(year_older) + '/' + str(month)
            if  (not str_year_month in acceptance_by_journal[journal_name].keys()):
                acceptance_by_journal[journal_name][str_year_month] = dict()
                acceptance_by_journal[journal_name][str_year_month]['total_days'] = 0
                acceptance_by_journal[journal_name][str_year_month]['num_articles'] = 0
                if journal_name.lower().endswith('the lancet'):
                    empty_results.add(journal_name)
                    continue
                res_journals = self.get_journal_data(journal_name, str_year_month, str_year_month)
                if int(res_journals['Count']) == 0:
                    tmp_journal_name = journal_name.replace("&", "")
                    res_journals = self.get_journal_data(tmp_journal_name, str_year_month, str_year_month)
                    # num_publication[journal_name][str_year_month]['total'] = str(res_journals['Count'])
                    if int(res_journals['Count']) == 0:
                        empty_results.add(journal_name + " " + str(month))

                print('working on {}'.format(journal_name))
                acceptance_by_journal = self.get_paper_dates(acceptance_by_journal, journal_name, res_journals,str_year_month,None)

        for k,v in acceptance_by_journal.items():
            print(k,v)
        print(empty_results)
        print(len(empty_results))
        return acceptance_by_journal
        # pss.extract_journal_data(res_journals, journal_name)

    def get_paper_dates(self, acceptance_by_journal, journal_name, res_journals, str_year_month, covid_acceptance_data):
        num_handled_articles = 0
        count_invalid = 0
        if covid_acceptance_data!=None:
            covid_ids=covid_acceptance_data['ids']
        for id in res_journals['IdList']:
            if covid_acceptance_data != None:
                if id in covid_ids:
                    continue
            res = self.fetch_data_from_pubmed(int(id))
            if not isinstance(res, str):
                if not isinstance(res['PubmedArticle'], list):
                    continue
                try:
                    num_handled_articles += 1

                    dates = res['PubmedArticle'][0]['PubmedData']['History']
                    for article_date in dates:
                        if article_date.attributes['PubStatus'] == 'received':
                            date_received = date(int(article_date['Year']), int(article_date['Month']),
                                                 int(article_date['Day']))
                        if article_date.attributes['PubStatus'] == 'revised':
                            date_revised = date(int(article_date['Year']), int(article_date['Month']),
                                                int(article_date['Day']))
                        if article_date.attributes['PubStatus'] == 'accepted':
                            date_accepted = date(int(article_date['Year']), int(article_date['Month']),
                                                 int(article_date['Day']))

                    acceptance_time = date_accepted - date_received + timedelta(days=1)
                    print(acceptance_time)
                    if acceptance_time.total_seconds() <= 0:
                        count_invalid += 1
                        continue
                    acceptance_by_journal[journal_name][str_year_month]['total_days'] += acceptance_time.days
                    acceptance_by_journal[journal_name][str_year_month]['num_articles'] += 1

                except:
                    print('caught exception for {}'.format(res))
        return acceptance_by_journal


if __name__ == '__main__':
    print(datetime.now())
    utils = Utils(path='D:\\shir\\study\\covid_19\\scopus')
    aa = AcceptanceAnalysis()
    vis=Visualization()
    # aa.iterate_journals()
    #
    # exit(0)
    data_dict=dict()
    months=['Jan']
    for month in months:
        data_dict[month]=dict()
        acc_speed, acceptance_by_journal, num_total_articles, num_handled_articles, num_invalid_speed_articles=aa.get_dates(month+'_articles.csv')
        print(acc_speed)
        data_dict[month]['Total_days']=num_handled_articles-num_invalid_speed_articles
        data_dict[month]['acc_speed']=acc_speed
        data_dict[month]['journals']=acceptance_by_journal
    aa.iterate_journals(data_dict)
    # vis.plt_acceptance_speed(data_dict,"acceptance speed")
    print(datetime.now())
