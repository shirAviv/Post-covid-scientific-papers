import os
import time
from bio import Entrez
from utils import Utils
import pickle
import pandas as pd
from bio import Entrez
import numpy as np
from datetime import date,datetime,timedelta


path='D:\\shir\\study\\covid_19\\scopus'
default_date = date(2030, 12, 12)


class SciencedirectData:

    def __init__(self):
        self.utils=Utils(path)

    def get_data(self,file):
        csv_path = os.path.join(path, file)
        papers = pd.read_csv(csv_path)
        papers = papers.dropna(how='any', subset=['DOI'])
        papers = papers.drop(['Journal'], axis=1)
        papers[['Date_Received', 'Date_Accepted', 'Date_online', 'Acceptance_Time', 'Countries', 'Pub_name', 'Document_type']] = papers.apply(
            lambda row: pd.Series(self.extract(row['DOI'], row['Title'])), axis=1)
        # contents = utils.read_from_csv(file).split('\n')
        # contents = contents[0:-1]
        # authors_dict = dict()
        # countries_totals_dict = dict()
        # countries_rel_dict = dict()
        # countries_per_article_dict = dict()
        # handled_current_article = False
        # unknowns = []
        # num_handled_articles = 0
        # for idx, line in enumerate(contents):
        #     if len(line) == 1:
        #         continue
        #     title, authors, doi, journal_name, date, url = line.split(',')
        # date_received,date_accepted, date_online, acceptance_time,countries, pub_name,pub_type = self.extract(doi, title)
        print("{} {}".format(papers['Pub_name'],papers['Date_online']))
        return papers


    def extract(self, doi, title):
        res, affils, pub_name, pub_type = self.utils.get_data_from_doi(doi,title)
        countries = set()
        if affils!=None:
            for affil in affils:
                countries.add(affil['affiliation-country'])
        date_received = np.nan
        date_accepted = np.nan
        acceptance_time = np.nan
        date_online=default_date
        if not isinstance(res, str):
            if isinstance(res['PubmedArticle'], list):
                date_received, date_accepted, date_online, acceptance_time = self.dates_calc(res)
            if pub_name == None:
                pub_name=self.extract_pub_name_from_pubmed(res)
        return date_received,date_accepted, date_online, acceptance_time,countries, pub_name, pub_type

    def extract_pub_name_from_pubmed(self,res):
        pub_name=None
        try:
            pub_name = res['PubmedArticle'][0]['MedlineCitation']['Article']['Journal']['Title']
        finally:
            return pub_name

    def dates_calc(self,res):

        date_received = default_date
        date_accepted = default_date
        date_online = default_date
        acceptance_time = date_accepted-date_received
        try:
            id = res['PubmedArticle'][0]['MedlineCitation']['PMID']
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
        finally:
            try:
                date_online = res['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleDate'][0]
                date_online = date(int(date_online['Year']), int(date_online['Month']), int(date_online['Day']))
            except:
                date_online = date_accepted
            return date_received, date_accepted, date_online, acceptance_time


if __name__ == '__main__':
    start_time = datetime.now()
    print(start_time)
    utils=Utils(path=path)
    pa=SciencedirectData()
    papers=pa.get_data('Apr_art1.csv')
    utils.write_to_csv(papers,os.path.join(path, 'Apr_articles_extended1' + '.csv'))
    end_time = datetime.now()
    print(end_time - start_time)
