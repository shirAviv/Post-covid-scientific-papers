import os
import time
from Bio import Entrez
from utils import Utils
import pickle
import pandas as pd
from Bio import Entrez
import numpy as np
from datetime import date,datetime,timedelta
from parse_scopus_search import ParseScopusSearch
from sciencedirect_data import SciencedirectData

path='D:\\shir\\study\\covid_19\\scopus'

class NumPublicationsHistory:

    def get_publications_for_year_month(self, stats, year,month=None):
        if month!=None:
            str_year_month=str(year)+'/'+str(month)
        else:
            str_year_month = str(year)
        stats[[str_year_month]] = stats.apply(
            lambda row: pd.Series(self.get_publications_num(row['Pub_name'],str_year_month)), axis=1)
        return stats

    def get_publications_num(self, pub_name,str_year_month):
            results = pss.get_journal_data(pub_name, str_year_month, str_year_month)
            return results['Count']

if __name__ == '__main__':
    pss=ParseScopusSearch()
    utils=Utils(path=path)
    nph=NumPublicationsHistory()
    journals=utils.load_csv_data_to_df('journals_history.csv')
    stats=journals
    for year in range(2016,2021):
        for month in range(1, 5):
            stats1 = nph.get_publications_for_year_month(journals, year, month)
            str_year_month = str(year) + '/' + str(month)
            stats[str_year_month] = stats1[str_year_month].values
        stats1=nph.get_publications_for_year_month(journals,year)
        stats[str(year)]=stats1[str(year)].values

    utils.write_to_csv(stats,'journals_history_1.csv')


    # res_journals = pss.iterate_journals()
    # res_journals = pss.handle_special_journals(res_journals)
    # utils.save_obj(res_journals,"journals_data_apr")
    exit(0)