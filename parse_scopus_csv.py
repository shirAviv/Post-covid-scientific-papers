import requests
import json
import pandas as pd
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


class ParseScopusCsv:

    def __init__(self):
        self.utils=Utils(path)

    def extract_covid_from_scopus(self, file):
        csv_path = os.path.join(path, file)
        papers=pd.read_csv(csv_path)
        papers=papers.dropna(how='any', subset=['DOI'])
        papers=papers.drop(['Volume','Issue','Art. No.','Page start','Page end', 'Page count'], axis=1)
        papers[['Date_Received','Date_Accepted','Date_online','Acceptance_Time','Countries']]=papers.apply(lambda row: pd.Series(self.get_dates_and_countries(row['DOI'],row['Title'])),axis=1)

        # papers=self.apply_and_concat(papers,['DOI','Title'], self.get_dates_and_countries,['Date_Received','Date_Accepted','Date_online','Acceptance_Time','Countries'])


        utils.save_obj(papers,"scopus_data_4000")

        # papers[['Date_Received','Date_Accepted','Date_online','Acceptance_Time','Countries']]=papers.apply(lambda row: pd.Series(self.get_dates_and_countries(row['DOI'])),axis=1)
        # papers_Feb=papers[papers['Date_online'].month==2]
        # for row in papers.itertuples(index=True, name='Pandas'):
        #     doi=getattr(row,'DOI')
        #     print(getattr(row, "Date_online"),getattr(row,'Countries'))

    def apply_and_concat(self,dataframe, field, func, column_names):
        return pd.concat((
            dataframe,
            dataframe[field].apply(
                lambda cell: pd.Series(func(cell), index=column_names))), axis=1)

    def get_dates_and_countries(self,doi,title):
        res, affils = self.utils.get_data_from_doi(doi,title)
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
        return date_received,date_accepted, date_online, acceptance_time,countries

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

    def join_and_save(self, papers1, papers2):
        papers=pd.concat([papers1,papers2],ignore_index=True)
        print(papers)
        print(len(papers))
        utils.save_obj(papers,'scopus_data')


    def split_by_month(self, papers):
        papers['Date_online'] = pd.to_datetime(papers['Date_online'])
        papers_Jan = papers[papers['Date_online'] < pd.Timestamp(2020, 2, 1)]
        print("papers in Jan {}".format(len(papers_Jan)))
        paper_Jan_articles=papers_Jan.loc[(papers_Jan['Document_Type'] == 'Article') | (papers_Jan['Document_Type'] == 'Review')]
        papers_Feb = papers[papers['Date_online'] < pd.Timestamp(2020, 3, 1)]
        papers_Feb = papers_Feb[papers_Feb['Date_online'] >= pd.Timestamp(2020, 2, 1)]
        print("papers in Feb {}".format(len(papers_Feb)))
        paper_Feb_articles=papers_Feb.loc[(papers_Feb['Document_Type'] == 'Article') | (papers_Feb['Document_Type'] == 'Review')]
        papers_Mar = papers[papers['Date_online'] < pd.Timestamp(2020, 4, 1)]
        papers_Mar = papers_Mar[papers_Mar['Date_online'] >= pd.Timestamp(2020, 3, 1)]
        print("papers in Mar {}".format(len(papers_Mar)))
        paper_Mar_articles=papers_Mar.loc[(papers_Mar['Document_Type'] == 'Article') | (papers_Mar['Document_Type'] == 'Review')]
        papers_Apr = papers[papers['Date_online'] < pd.Timestamp(2020, 5, 1)]
        papers_Apr = papers_Apr[papers_Apr['Date_online'] >= pd.Timestamp(2020, 4, 1)]
        print("papers in Apr {}".format(len(papers_Apr)))
        paper_Apr_articles=papers_Apr.loc[(papers_Apr['Document_Type'] == 'Article') | (papers_Apr['Document_Type'] == 'Review')]
        papers_May = papers[papers['Date_online'] < pd.Timestamp(2020, 6, 1)]
        papers_May = papers_May[papers_May['Date_online'] >= pd.Timestamp(2020, 5, 1)]
        print("papers in May {}".format(len(papers_May)))
        paper_May_articles=papers_May.loc[(papers_May['Document_Type'] == 'Article') | (papers_May['Document_Type'] == 'Review')]
        print('articles and reviews Jan {}, Feb {}, Mar {}, Apr {}, May {}'.format(len(paper_Jan_articles),len(paper_Feb_articles),len(paper_Mar_articles), len(paper_Apr_articles), len(paper_May_articles)))
        return papers_Jan, papers_Feb, papers_Mar, papers_Apr, papers_May

    def get_publication_data(self,month,papers):
        sorted_papers=papers.sort_values(by='Source_title')
        # papers['publication_freq']=papers.groupby('Source_title')['Source_title'].transform('count')
        publications_counts=papers['Source_title'].value_counts()
        print(publications_counts)




if __name__ == '__main__':
    start_time=datetime.now()
    print(start_time)
    psc=ParseScopusCsv()
    utils=Utils(path=path)
    # a=utils.load_obj("scopus_data_6000")
    # b=utils.load_obj("scopus_data_4000")
    papers=utils.load_obj('scopus_data')
    pj,pf,pm,pa,pmy=psc.split_by_month(papers)
    psc.get_publication_data('jan',pj)
    psc.get_publication_data('feb',pf)

    # psc.join_and_save(a,b)
    exit(0)
    file='scopus_tmp.csv'
    psc.extract_covid_from_scopus(file)
    end_time=datetime.now()
    print(end_time-start_time)

