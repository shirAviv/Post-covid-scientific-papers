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
# from sciencedirect_data import SciencedirectData


path='D:\\shir\\study\\covid_19\\scopus'


class AcceptanceTime:
    def get_publication_data(self,month,papers):

        papers[['Acceptance_Time']] = papers.apply(
            lambda row: pd.Series(self.get_acceptance_time(row['Acceptance_Time'])), axis=1)

        # sorted_papers=papers.sort_values(by='Pub_name')
        group_by_pub_name=papers.groupby(['Pub_name'])

        stats=pd.DataFrame(columns=['Total_counts','count_nonzero','mean','std'])
        stats['Total_counts']=group_by_pub_name['Acceptance_Time'].agg(np.count_nonzero)
        stats['count_nonzero'] = group_by_pub_name['Acceptance_Time'].count()
        stats['mean']= group_by_pub_name['Acceptance_Time'].mean()
        stats['std']= group_by_pub_name['Acceptance_Time'].std()

        # stats1=group_by_pub_name['Acceptance_Time'].agg([np.count_nonzero, np.mean,np.std,np.median])
        stats.reset_index(inplace=True)
        stats=stats.sort_values(by='count_nonzero', ascending=False)[0:30]
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        print('stats for mon {} are {}'.format(month,stats))
        return stats

    def extract_and_apply_publicattion_metrics(self,file,stats):
        metrics=self.load_data(file)
        metrics['SJR']=pd.to_numeric(metrics['SJR'],errors='coerce')
        metrics['CiteScore']=pd.to_numeric(metrics['CiteScore'],errors='coerce')
        metrics['SNIP']=pd.to_numeric(metrics['SNIP'],errors='coerce')

        metrics=metrics.sort_values(by='SJR')
        stats[['SJR','CiteScore','SNIP']] = stats.apply(
            lambda row: pd.Series(self.get_metrics(row['Pub_name'],metrics)), axis=1)
        # for pub_name in stats.index.values:
        #     print(pub_name)
            # if not pub_name in metrics.values:
            #     print("{} not exists".format(pub_name))
            # journal_metrics=metrics[metrics['journal_name']==pub_name]
        print(stats)
        return stats

    def get_metrics(self,journal,metrics):
        try:
            data=metrics[metrics['journal_name']==journal]
            SJR=float(data['SJR'])
            CiteScore = float(data['CiteScore'])
            SNIP = float(data['SNIP'])
        except:
            SJR = 0
            CiteScore = 0
            SNIP = 0
        return SJR,CiteScore,SNIP


    def get_acceptance_time(self,accept_time):
        if accept_time!='NaT' and accept_time!=np.nan and accept_time!='':
            days= int(accept_time.split(' ')[0])
            if days<0:
                return np.nan
            else:
                return days
        else:
            return np.nan

    def load_data(self,file):
        csv_path = os.path.join(path, file)
        papers = pd.read_csv(csv_path,keep_default_na=False)
        return papers

    # def remove_dupes(self,papers):
    #     print(len(papers))
    #     papers1=papers.drop_duplicates(subset='DOI')
    #     print(len(papers1))

    def get_journals_by_month(self, covid_stats,month,year):
        str_year_month = str(year) + '/' + str(month)
        journals=covid_stats['Pub_name']
        for journal in journals:
            results = pss.get_journal_data(journal, str_year_month,str_year_month)
            print(results)
            df=pd.DataFrame(columns=['Pub_name','date_received','date_accepted', 'date_online', 'Acceptance_Time'])
            if int(results['Count']) > 0:
                ids = results['IdList']
            else:
                ids=None
            if ids != None:
                for index,id in enumerate(ids):
                    # if index==100:
                    #     break
                    paper=utils.fetch_data_from_pubmed(id)
                    # print(paper)
                    if not isinstance(paper, str):
                        if isinstance(paper['PubmedArticle'], list):
                            date_received,date_accepted,date_online,acceptance_time = scd.dates_calc(paper)
                            df=df.append({'Pub_name':journal,'date_received': date_received,'date_accepted': date_accepted,'date_online': date_online,'Acceptance_Time':acceptance_time}, ignore_index=True)
                df = df.astype({"Acceptance_Time": str})
                stats=self.get_publication_data(month,df)
                covid_stats.loc[covid_stats['Pub_name'] == journal, 'count_nonzero_'+str(year)+str(month)] = stats['count_nonzero'][0]
                covid_stats.loc[covid_stats['Pub_name'] == journal, 'mean_'+str(year)+str(month)] = stats['mean'][0]
                covid_stats.loc[covid_stats['Pub_name'] == journal, 'std_'+str(year)+str(month)] = stats['std'][0]
                covid_stats.loc[covid_stats['Pub_name'] == journal, 'median_'+str(year)+str(month)] = stats['median'][0]
        # print(covid_stats)
        return(covid_stats)


    def extract_stats_by_month(self,month,month_str):
        papers = test.load_data(month_str + '_articles_acceptance.csv')
        # test.remove_dupes(papers)
        stats = self.get_publication_data(month, papers)
        stats = self.extract_and_apply_publicattion_metrics('journals_list_metrics.csv', stats)
        curr_month_papers=stats
        # curr_month_papers = test.get_journals_by_month(stats, month, 2019)
        return curr_month_papers


if __name__ == '__main__':
    start_time = datetime.now()
    print(start_time)
    utils = Utils(path=path)
    test = AcceptanceTime()
    pss = ParseScopusSearch()
    # scd = SciencedirectData()
    full_stats=pd.DataFrame()
    # month=1
    # month_str='Jan'
    # curr_month_papers1=test.extract_stats_by_month(month,month_str)
    # month = 2
    # month_str = 'Feb'
    # curr_month_papers2 = test.extract_stats_by_month(month, month_str)
    # month = 3
    # month_str = 'Mar'
    # curr_month_papers3 = test.extract_stats_by_month(month, month_str)
    # month = 4
    # month_str = 'Apr'
    # curr_month_papers4 = test.extract_stats_by_month(month, month_str)
    # full_stats = pd.concat([curr_month_papers1, curr_month_papers2, curr_month_papers3,curr_month_papers4], ignore_index=True)
    # utils.write_to_csv(full_stats, 'covid_acceptance_time.csv')

    # stats=pd.DataFrame()
    year=2020
    # month = 1
    # month_str = 'Jan'
    # papers = test.load_data(month_str + '_articles_extended.csv')
    # stats['Pub_name']=papers['Pub_name'].unique()
    # curr_month_papers1 = test.get_journals_by_month(stats, month, 2020)
    # month = 2
    # month_str = 'Feb'
    # papers = test.load_data(month_str + '_articles_extended.csv')
    # stats = stats.iloc[0:0]
    # stats['Pub_name'] = papers['Pub_name'].unique().copy()
    # curr_month_papers2 = test.get_journals_by_month(stats, month, 2020)
    # month = 3
    # month_str = 'Mar'
    # papers = test.load_data(month_str + '_articles_extended.csv')
    # stats = stats.iloc[0:0]
    # stats['Pub_name'] = papers['Pub_name'].unique().copy()
    # curr_month_papers3 = test.get_journals_by_month(stats, month, 2020)

    month = 1
    month_str = 'Jan'
    stats=utils.load_csv_data_to_df('journals_history.csv')

    stats = test.get_journals_by_month(stats, month, 2020)
    # month=2
    # stats = test.get_journals_by_month(stats, month, 2020)
    # month = 3
    # stats = test.get_journals_by_month(stats, month, 2020)
    # month = 4
    # stats = test.get_journals_by_month(stats, month, 2020)


    # utils.write_to_csv(stats, '2020_acceptance_time.csv')
    print(full_stats)