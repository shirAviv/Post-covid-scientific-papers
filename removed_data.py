import os
import time
from bio import Entrez
from utils import Utils
import pickle
import pandas as pd
from bio import Entrez
import numpy as np
from sciencedirect_data import SciencedirectData
from datetime import date,datetime,timedelta
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import convert_to_unicode
from AcceptanceTime import AcceptanceTime





path='D:\\shir\\study\\covid_19\\scopus\scienceDirectData'

class RemovedData:

    def iterate_bibs(self):
        for file in os.listdir(path=path):
            if file.startswith("ScienceDirect_citations_") and file.endswith(".bib"):
                name=file.split('citations_')[1].split('.bib')[0]
                print(name)
                obj_path = os.path.join(path, file)
                bibtex_file = open(obj_path, "r", encoding="utf8")
                # read content of file to string
                data = bibtex_file.read()
                # get number of occurrences of the substring in the string
                occurrences = data.count("@art")
                print(occurrences)

                csv_path=os.path.join(path, name+".csv")
                df = utils.load_csv_data_to_df(csv_path)

                print(len(df))

                # papers = sdd.get_data(obj_path)
                # print(name,papers)



    def get_journals_counts_by_year_and_month(self, utils):
        df_counts = pd.DataFrame()
        counts_removed_dict=dict()
        for file in os.listdir(path=path):
            if file.endswith("_history.csv"):
                journal_name = file.replace('_history.csv', '')
                journal_name = journal_name.replace('_', ' ')
                file_covid = file.replace('_history', '_COVID')
                if journal_name=='New Scientist':
                    continue
                row_counts = dict()
                row_counts['pub_name'] = [journal_name + '_counts']

                row_accept_mean = dict()
                row_accept_mean['pub_name'] = [journal_name + '_acc_mean']

                df_history = utils.load_csv_data_to_df(file)
                df_history.drop_duplicates(subset=['DOI'], inplace=True)
                df_COVID = utils.load_csv_data_to_df(file_covid)
                df_COVID.drop_duplicates(subset=['DOI'], inplace=True)

                all_DOIS = list(df_history['DOI'])
                vals = df_COVID['DOI'].isin(all_DOIS)
                df_COVID = df_COVID.loc[vals == True]
                removed_journal=0

                for year in range(2016, 2021):
                    removed_year=0
                    total_year=0
                    str_year = str(year)
                    df_year = df_history[df_history['Date_online'].str.match('^' + str_year + '.*') == True]
                    df_year_COVID = df_COVID[df_COVID['Date_online'].str.match('^' + str_year + '.*') == True]

                    if not str_year in counts_removed_dict.keys():
                        counts_removed_dict[str_year]=dict()
                        counts_removed_dict[str_year]['_removed']=0
                        counts_removed_dict[str_year]['_total']=0
                    for month in range(1, 7):
                        str_month = '0' + str(month)
                        df_month = df_year[df_year['Date_online'].str.match('.*-' + str_month + '-.*') == True].copy()
                        df_month['Acceptance_Time'] = df_month.apply(
                            lambda row: pd.Series(at.get_acceptance_time(row['Acceptance_Time'])), axis=1)

                        df_month_COVID = df_year_COVID[
                            df_year_COVID['Date_online'].str.match('.*-' + str_month + '-.*') == True].copy()
                        df_month_COVID['Acceptance_Time'] = df_month_COVID.apply(
                            lambda row: pd.Series(at.get_acceptance_time(row['Acceptance_Time'])), axis=1)

                        COVID_dois = list(df_month_COVID['DOI'])
                        vals = df_month['DOI'].isin(COVID_dois)
                        df_month_Non_COVID = df_month.loc[vals == False]

                        row_counts[str_year + '-' + str_month] = [df_month['DOI'].count()]
                        row_accept_mean[str_year + '-' + str_month] = [df_month['Acceptance_Time'].mean()]
                        if (len(df_month) > 0):
                            num_valid=len(df_month.loc[np.isnan(df_month['Acceptance_Time']) == False])
                            total=df_month['DOI'].count()
                            percent=100*num_valid/total
                        else:
                            num_valid=0
                            total=0
                            percent=100
                        # print(journal_name,year,month)
                        counts_removed = total - num_valid
                        total_year+=total
                        # print(counts_removed)
                        removed_year+=counts_removed
                        removed_journal+=counts_removed

                        # print(total)
                        # print(num_valid)
                        # print('{}%'.format(percent))
                        # if percent<80:
                        #     print("percent lower then 80")


                        if year == 2020:
                            row_counts[str_year + '-' + str_month + '_COVID'] = [df_month_COVID['DOI'].count()]
                            row_counts[str_year + '-' + str_month + '_Non_COVID'] = [df_month_Non_COVID['DOI'].count()]

                            row_accept_mean[str_year + '-' + str_month + '_COVID'] = [
                                df_month_COVID['Acceptance_Time'].mean()]
                            row_accept_mean[str_year + '-' + str_month + '_Non_COVID'] = [
                                df_month_Non_COVID['Acceptance_Time'].mean()]

                            if (len(df_month_COVID)>0):
                                num_valid = len(df_month_COVID.loc[np.isnan(df_month_COVID['Acceptance_Time']) == False])
                                total = df_month_COVID['DOI'].count()
                                percent = 100 * num_valid / total
                            else:
                                num_valid=0
                                total = 0
                                percent = 100
                            # print('covid')
                            # print(total-num_valid)

                            # print(total)
                            # print(num_valid)
                            # print('{}%'.format(percent))
                            # if percent < 80:
                            #     print("percent lower then 80")

                            if (len(df_month_Non_COVID) > 0):
                                num_valid = len(df_month_Non_COVID.loc[np.isnan(df_month_Non_COVID['Acceptance_Time']) == False])
                                total = df_month_Non_COVID['DOI'].count()
                                percent = 100 * num_valid / total
                            else:
                                num_valid=0
                                total = 0
                                percent = 100
                            # print('non covid')
                            # print(total-num_valid)
                            # print(total)
                            # print(num_valid)
                            # print('{}%'.format(percent))
                            # if percent < 80:
                            #     print("percent lower then 80")
                    counts_removed_dict[str_year]['_removed']+=removed_year
                    counts_removed_dict[str_year]['_total']+=total_year
                if not journal_name in counts_removed_dict.keys():
                    counts_removed_dict[journal_name]=0
                counts_removed_dict[journal_name]=removed_journal
        for k,v in counts_removed_dict.items():
            print(k,v)
            if k.startswith('20'):
                per=100*v['_removed']/v['_total']
                print('Percent {}%'.format(per))




if __name__ == '__main__':
    start_time = datetime.now()
    print(start_time)
    utils=Utils(path=path)
    sdd=SciencedirectData()
    at=AcceptanceTime()
    t=RemovedData()
    # t.iterate_bibs()
    t.get_journals_counts_by_year_and_month(utils)
