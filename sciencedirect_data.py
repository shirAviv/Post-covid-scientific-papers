import os
import time
from bio import Entrez
from utils import Utils
import pickle
import pandas as pd
from bio import Entrez
import numpy as np
from datetime import date,datetime,timedelta
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import convert_to_unicode
from AuthorsAndCountries import AuthorsAndCountries
from AcceptanceTime import AcceptanceTime
from visualization import Visualization




path='D:\\shir\\study\\covid_19\\scopus\scienceDirectData'
default_date = date(2030, 12, 12)


class SciencedirectData:

    def __init__(self):
        self.utils=Utils(path)

    def get_data(self,file):
        # csv_path = os.path.join(path, file)
        # papers = pd.read_csv(csv_path)
        papers=self.parse_bib_get_DOIs(file)
        papers = papers.dropna(how='any', subset=['DOI'])
        # papers = papers.drop(['Journal'], axis=1)
        papers[['Date_Received', 'Date_Accepted', 'Date_online', 'Acceptance_Time', 'Countries']] = papers.apply(
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
        print("{} {}".format(papers['journal'],papers['Date_online']))
        return papers

    def get_countries_from_res(self,results, doi):
        authors_dict = dict()
        countries_totals_dict = dict()
        countries_rel_dict = dict()
        countries_by_num_articles_dict = dict()
        unknowns = []
        num_handled_articles = 0
        articles_count = 0
        try:
            if len(results['PubmedArticle']) > 0 and ('MedlineCitation' in results['PubmedArticle'][0].keys()) and (
                    'Article' in results['PubmedArticle'][0]['MedlineCitation'].keys()):
                if 'AuthorList' in results['PubmedArticle'][0]['MedlineCitation']['Article'].keys():
                    authors_list = results['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']
                else:
                    print("no authors list {}".format(results['PubmedArticle'][0]['MedlineCitation']['Article']))
                    return None
                countries_by_num_articles_dict, countries_for_article, countries_totals_dict, \
                num_handled_articles, unknowns = aaa.get_authors_and_countries(
                    authors_dict, countries_by_num_articles_dict, countries_rel_dict, countries_totals_dict, doi,
                    num_handled_articles, authors_list, unknowns)
                return countries_for_article
        except:
            print("results are {}".format(results))
            return None


    def extract(self, doi, title):
        res, affils, pub_name, pub_type = self.utils.get_data_from_doi(doi,title)
        countries = set()
        if affils!=None:
            for affil in affils:
                countries.add(affil['affiliation-country'])
        else:
            countries=self.get_countries_from_res(res, doi)
        date_received = np.nan
        date_accepted = np.nan
        acceptance_time = np.nan
        date_online=default_date
        if not isinstance(res, str):
            if isinstance(res['PubmedArticle'], list):
                date_received, date_accepted, date_online, acceptance_time = self.dates_calc(res)
            # if pub_name == None:
            #     pub_name=self.extract_pub_name_from_pubmed(res)
        return date_received,date_accepted, date_online, acceptance_time,countries

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
            # acceptance_time = date_accepted - date_received + timedelta(days=1)
        finally:
            try:
                date_online = res['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleDate'][0]
                date_online = date(int(date_online['Year']), int(date_online['Month']), int(date_online['Day']))
            except:
                date_online = date_accepted
            if date_accepted==default_date and date_received!=default_date:
                date_accepted=date_online
            acceptance_time = date_accepted - date_received + timedelta(days=1)
            if acceptance_time.days<0:
                acceptance_time=np.nan
            if date_received==default_date:
                acceptance_time=np.nan
            return date_received, date_accepted, date_online, acceptance_time

    def join(self, papers1, papers2):
        papers = pd.concat([papers1, papers2], ignore_index=True)
        print(papers)
        print(len(papers))
        return papers

    def parse_bib_get_DOIs(self, file):
        with open(file, encoding="utf8") as bibtex_file:
            parser = BibTexParser()
            parser.customization = convert_to_unicode
            bib_database = bibtexparser.load(bibtex_file, parser=parser)
            print(bib_database.entries)
        df=pd.DataFrame(columns=['DOI','Title','journal','doc_type','url','authors'])
        for idx,paper in enumerate(bib_database.entries):
            doi=paper['doi']
            if not 'title' in paper.keys():
                print("missing title from doi {}".format(doi))
                continue
            Title=paper['title']
            journal = paper['journal']
            doc_type=paper['ENTRYTYPE']
            url=paper['url']
            if('author') in paper.keys():
                authors=paper['author']
            else:
                authors=None
            df.loc[idx]=[doi, Title, journal, doc_type, url, authors]
        return(df)

    def test_covid_in_data(self, history_df, covid_df):
        print('size of covid - {}'.format(len(covid_df)))
        print('size of all - {}'.format(len(history_df)))

        for row in covid_df.iterrows():
            doi=row[1]['DOI']
            val=len(history_df.loc[history_df['DOI']==doi])
            if val==0:
                print('doi not found {}'.format(doi))

    def get_journals_counts_by_year_and_month(self):

        df_counts=pd.DataFrame()
        for file in os.listdir(path=path):
            if file.endswith("_history.csv"):
                journal_name=file.replace('_history.csv','')
                journal_name=journal_name.replace('_',' ')
                file_covid=file.replace('_history','_COVID')

                row_counts=dict()
                row_counts['pub_name']=[journal_name+'_counts']


                row_accept_mean = dict()
                row_accept_mean['pub_name'] = [journal_name + '_acc_mean']

                df_history = utils.load_csv_data_to_df(file)
                df_history.drop_duplicates(subset=['DOI'], inplace=True)
                df_COVID =utils.load_csv_data_to_df(file_covid)
                df_COVID.drop_duplicates(subset=['DOI'], inplace=True)

                all_DOIS = list(df_history['DOI'])
                vals = df_COVID['DOI'].isin(all_DOIS)
                df_COVID = df_COVID.loc[vals == True]

                for year in range (2016,2021):
                    str_year = str(year)
                    df_year=df_history[df_history['Date_online'].str.match('^'+str_year+'.*')==True]
                    df_year_COVID=df_COVID[df_COVID['Date_online'].str.match('^'+str_year+'.*')==True]

                    for month in range(1,7):
                        str_month = '0' + str(month)
                        df_month = df_year[df_year['Date_online'].str.match('.*-' + str_month + '-.*') == True].copy()
                        df_month['Acceptance_Time'] = df_month.apply(
                            lambda row: pd.Series(at.get_acceptance_time(row['Acceptance_Time'])), axis=1)

                        df_month_COVID = df_year_COVID[df_year_COVID['Date_online'].str.match('.*-' + str_month + '-.*') == True].copy()
                        df_month_COVID['Acceptance_Time'] = df_month_COVID.apply(
                            lambda row: pd.Series(at.get_acceptance_time(row['Acceptance_Time'])), axis=1)

                        COVID_dois=list(df_month_COVID['DOI'])
                        vals=df_month['DOI'].isin(COVID_dois)
                        df_month_Non_COVID= df_month.loc[vals == False]

                        row_counts[str_year+'-'+str_month]=[df_month['DOI'].count()]
                        row_accept_mean[str_year+'-'+str_month]=[df_month['Acceptance_Time'].mean()]

                        if year==2020:
                            row_counts[str_year+'-'+str_month+'_COVID']=[df_month_COVID['DOI'].count()]
                            row_counts[str_year + '-' + str_month+'_Non_COVID'] = [df_month_Non_COVID['DOI'].count()]

                            row_accept_mean[str_year+'-'+str_month+'_COVID']=[df_month_COVID['Acceptance_Time'].mean()]
                            row_accept_mean[str_year+'-'+str_month+'_Non_COVID']=[df_month_Non_COVID['Acceptance_Time'].mean()]



                df_counts=df_counts.append(pd.DataFrame.from_dict(row_counts))
                df_counts=df_counts.append(pd.DataFrame.from_dict(row_accept_mean))


        print(df_counts)

        return df_counts


    def extract_countries_data(self):
        # rank countries by num publications, each year, each month
        # rank countries by num COVID publications, each month
        #  internationally collaborative papers - in a given month and year, for all ccountries :count number of papers with one author, with cross country two, three.. authors.
            # compare to COVID
        # collaboration - avg # of authors per paper in each country, per month per year and compared to avg
            # of authors per COVID paper in each country

        df_countries = pd.DataFrame()
        countries_dict=dict()
        for file in os.listdir(path=path):
            if file.endswith("_history.pkl"):
                file_to_load=file.replace('.pkl','')
                journal_name = file.replace('_history.pkl', '')
                journal_name = journal_name.replace('_', ' ')
                file_covid = file_to_load.replace('_history', '_COVID')

                row_counts = dict()
                row_counts['pub_name'] = [journal_name + '_counts']


                df_history = utils.load_obj(file_to_load)
                df_history.drop_duplicates(subset=['DOI'], inplace=True)
                df_history=df_history.dropna()
                df_history = df_history.astype({'Date_online': str})

                df_COVID = utils.load_obj(file_covid)
                df_COVID.drop_duplicates(subset=['DOI'], inplace=True)
                df_COVID=df_COVID.dropna()
                df_COVID = df_COVID.astype({'Date_online': str})

                for year in range (2016,2021):
                    if not str(year) in countries_dict.keys():
                        countries_dict[str(year)]=dict()


                    str_year = str(year)
                    df_year=df_history[df_history['Date_online'].str.match('^'+str_year+'.*')==True]
                    df_year_COVID=df_COVID[df_COVID['Date_online'].str.match('^'+str_year+'.*')==True]

                    for month in range(1,7):
                        if not str(month) in countries_dict[str(year)].keys():
                            countries_dict[str(year)][str(month)] = dict()
                        str_month = '0' + str(month)
                        df_month = df_year[df_year['Date_online'].str.match('.*-' + str_month + '-.*') == True].copy()

                        df_month_COVID = df_year_COVID[df_year_COVID['Date_online'].str.match('.*-' + str_month + '-.*') == True].copy()

                        COVID_dois=list(df_month_COVID['DOI'])
                        vals=df_month['DOI'].isin(COVID_dois)
                        df_month_Non_COVID= df_month.loc[vals == False]
                        if year==2020:
                            if not str(month)+'_COVID' in countries_dict[str(year)].keys():
                                countries_dict[str(year)][str(month)+'_COVID'] = dict()

                        countries=df_month.Countries
                        if not 'Monthly Total articles' in countries_dict[str(year)][str(month)]:
                            countries_dict[str(year)][str(month)]['Monthly Total articles'] = 0
                        countries_dict[str(year)][str(month)]['Monthly Total articles']+=len(countries)
                        for country_set in countries:
                            for country in country_set:
                                if country=='unknown':
                                    continue
                                try:
                                    subset_without_country=country_set.copy()
                                    subset_without_country.remove(country)
                                except:
                                    print(country_set)
                                    continue
                                if not country in countries_dict[str(year)][str(month)].keys():
                                    countries_dict[str(year)][str(month)][country]=dict()
                                    countries_dict[str(year)][str(month)][country]['Country Total Articles']=0
                                countries_dict[str(year)][str(month)][country]['Country Total Articles'] += 1
                                for curr_country in subset_without_country:
                                    if curr_country == 'unknown':
                                        continue
                                    if not curr_country in countries_dict[str(year)][str(month)][country].keys():
                                        countries_dict[str(year)][str(month)][country][curr_country]=0
                                    countries_dict[str(year)][str(month)][country][curr_country] += 1

                        if year==2020:
                            countries = df_month_COVID.Countries
                            if not 'Monthly Total articles' in countries_dict[str(year)][str(month)+'_COVID']:
                                countries_dict[str(year)][str(month) + '_COVID']['Monthly Total articles'] = 0
                            countries_dict[str(year)][str(month)+'_COVID']['Monthly Total articles'] += len(countries)
                            for country_set in countries:
                                for country in country_set:
                                    if country == 'unknown':
                                        continue
                                    try:
                                        subset_without_country = country_set.copy()
                                        subset_without_country.remove(country)
                                    except:
                                        print(country_set)
                                        continue
                                    if not country in countries_dict[str(year)][str(month)+'_COVID'].keys():
                                        countries_dict[str(year)][str(month)+'_COVID'][country] = dict()
                                        countries_dict[str(year)][str(month)+'_COVID'][country]['Country Total Articles'] = 0
                                    countries_dict[str(year)][str(month)+'_COVID'][country]['Country Total Articles'] += 1
                                    for curr_country in subset_without_country:
                                        if curr_country == 'unknown':
                                            continue
                                        if not curr_country in countries_dict[str(year)][str(month)+'_COVID'][country].keys():
                                            countries_dict[str(year)][str(month)+'_COVID'][country][curr_country] = 0
                                        countries_dict[str(year)][str(month)+'_COVID'][country][curr_country] += 1

                        # countries_year_dict[str(year)][str(month)]=countries_month


                print(journal_name)
        return countries_dict

    def get_list_from_set(self, countries):
        return list(countries)

    def create_counts_table(self,df):
        df_counts=df[df['pub_name'].str.match('.*_counts')].copy()
        df_counts.set_index('pub_name', inplace=True)
        df_counts.rename(index={0:'pub_name'}, inplace=True)
        df_counts = df_counts.apply(pd.to_numeric, errors='coerce')
        # df_counts['SJR']=[]
        df_counts_COVID=df_counts.loc[:, '2020-01':'2020-06_Non_COVID'].copy()
        df_counts_COVID.reset_index(inplace=True)
        df_counts_COVID.rename(columns={'index': 'pub_name'}, inplace=True)
        metrics = utils.load_csv_data_to_df('journals_list_metrics_tmp.csv')
        metrics['SJR'] = pd.to_numeric(metrics['SJR'], errors='coerce')
        metrics['CiteScore'] = pd.to_numeric(metrics['CiteScore'], errors='coerce')
        metrics['SNIP'] = pd.to_numeric(metrics['SNIP'], errors='coerce')
        df_counts_COVID['pub_name'] = df_counts_COVID.apply(
            lambda row: pd.Series(row['pub_name'].replace('_counts', '')), axis=1)
        df_counts_COVID['SJR'] = df_counts_COVID.apply(
            lambda row: pd.Series(self.get_SJR(row['pub_name'], metrics)), axis=1)


        df_counts_COVID['Percent_01'] = df_counts_COVID.apply(lambda row: '{:.2%}'.format(row['2020-01_COVID'] / row['2020-01']), axis=1)
        df_counts_COVID['Percent_02'] = df_counts_COVID.apply(lambda row: '{:.2%}'.format(row['2020-02_COVID'] / row['2020-02']), axis=1)
        df_counts_COVID['Percent_03'] = df_counts_COVID.apply(lambda row: '{:.2%}'.format(row['2020-03_COVID'] / row['2020-03']), axis=1)
        df_counts_COVID['Percent_04'] = df_counts_COVID.apply(lambda row: '{:.2%}'.format(row['2020-04_COVID'] / row['2020-04']), axis=1)
        df_counts_COVID['Percent_05'] = df_counts_COVID.apply(lambda row: '{:.2%}'.format(row['2020-05_COVID'] / row['2020-05']), axis=1)
        df_counts_COVID['Percent_06'] = df_counts_COVID.apply(lambda row: '{:.2%}'.format(row['2020-06_COVID'] / row['2020-06']), axis=1)

        print(df_counts_COVID)
        df_counts_COVID=df_counts_COVID.reindex(columns=['pub_name','SJR','2020-01_COVID','2020-01',
                                             '2020-02_COVID','2020-02',
                                             '2020-03_COVID','2020-03',
                                             '2020-04_COVID','2020-04',
                                             '2020-05_COVID','2020-05',
                                             '2020-06_COVID','2020-06',
                                            'Percent_01', 'Percent_02',
                                            'Percent_03', 'Percent_04',
                                            'Percent_05', 'Percent_06',])
        df_counts_COVID = df_counts_COVID.sort_values('SJR', ascending=False, na_position='last')
        df_counts_COVID.set_index('pub_name', inplace=True)
        # df_counts_COVID=df_counts_COVID.loc[:,[]]
        # df_counts_COVID.style.set_table_styles([{'selector': '', 'props': [('border', '4px solid #7a7')]}])
        tex=df_counts_COVID.to_latex()
        text_file = open("covid_growth_table.tex", "w")
        text_file.write(tex)
        text_file.close()

    def get_history_counts(self, df):
        df_counts=df[df['pub_name'].str.match('.*_counts')].copy()
        df_counts.set_index('pub_name', inplace=True)
        df_counts.rename(index={0:'pub_name'}, inplace=True)
        df_counts=df_counts.drop(['New Scientist_counts'])
        df_counts = df_counts.apply(pd.to_numeric, errors='coerce')
        all_counts_sum=pd.DataFrame()
        all_counts_sum[0] = df_counts.sum()
        all_counts_sums_history=all_counts_sum.loc['2016-01':'2020-01',:]
        all_counts_sums_history=all_counts_sums_history.append(all_counts_sum.loc[['2020-02', '2020-03','2020-04','2020-05','2020-06'],:])
        all_counts_sums_history.rename(columns={0:'Total'}, inplace=True)
        covId_counts_sums=pd.DataFrame(index=all_counts_sums_history.index)
        covId_counts_sums.loc['2020-01',0]=all_counts_sum.loc['2020-01_COVID',0]
        covId_counts_sums.loc['2020-02',0] = all_counts_sum.loc['2020-02_COVID', 0]
        covId_counts_sums.loc['2020-03',0] = all_counts_sum.loc['2020-03_COVID', 0]
        covId_counts_sums.loc['2020-04',0] = all_counts_sum.loc['2020-04_COVID', 0]
        covId_counts_sums.loc['2020-05',0] = all_counts_sum.loc['2020-05_COVID', 0]
        covId_counts_sums.loc['2020-06',0] = all_counts_sum.loc['2020-06_COVID', 0]
        # all_counts_sum.loc[['2020-01_COVID','2020-02_COVID', '2020-03_COVID','2020-04_COVID','2020-05_COVID','2020-06_COVID'],:]
        covId_counts_sums.rename(columns={0:'Total'}, inplace=True)
        vis.plt_journals_sums(all_counts_sums_history, covId_counts_sums,  "COVID-19 and Total publications, Scopus journals")


    def get_acc_time(self,df):
        df_means = df[df['pub_name'].str.match('.*_mean')].copy()
        df_means.set_index('pub_name', inplace=True)
        df_means.rename(index={0: 'pub_name'}, inplace=True)
        df_means = df_means.drop(['New Scientist_acc_mean'])
        df_means = df_means.apply(pd.to_numeric, errors='coerce')
        df_means_avg=pd.DataFrame()
        df_means_avg[0] = df_means.mean()
        df_means_sem=pd.DataFrame()
        df_means_sem[0]= df_means.sem()
        df_means_avg_history = df_means_avg.loc['2016-01':'2020-01']
        df_means_avg_history = df_means_avg_history.append(
            df_means_avg.loc[['2020-02', '2020-03', '2020-04', '2020-05', '2020-06']])
        df_means_sem_history = df_means_sem.loc['2016-01':'2020-01']
        df_means_sem_history = df_means_sem_history.append(
            df_means_sem.loc[['2020-02', '2020-03', '2020-04', '2020-05', '2020-06']])

        df_means_avg_history.rename(columns={0: 'Total'}, inplace=True)
        df_means_sem_history.rename(columns={0: 'Total'}, inplace=True)

        covId_means_avg = pd.DataFrame(index=df_means_avg_history.index)
        covId_means_sem = pd.DataFrame(index=df_means_sem_history.index)

        covId_means_avg.loc['2020-01', 0] = df_means_avg.loc['2020-01_COVID', 0]
        covId_means_avg.loc['2020-02', 0] = df_means_avg.loc['2020-02_COVID', 0]
        covId_means_avg.loc['2020-03', 0] = df_means_avg.loc['2020-03_COVID', 0]
        covId_means_avg.loc['2020-04', 0] = df_means_avg.loc['2020-04_COVID', 0]
        covId_means_avg.loc['2020-05', 0] = df_means_avg.loc['2020-05_COVID', 0]
        covId_means_avg.loc['2020-06', 0] = df_means_avg.loc['2020-06_COVID', 0]

        covId_means_sem.loc['2020-01', 0] = df_means_sem.loc['2020-01_COVID', 0]
        covId_means_sem.loc['2020-02', 0] = df_means_sem.loc['2020-02_COVID', 0]
        covId_means_sem.loc['2020-03', 0] = df_means_sem.loc['2020-03_COVID', 0]
        covId_means_sem.loc['2020-04', 0] = df_means_sem.loc['2020-04_COVID', 0]
        covId_means_sem.loc['2020-05', 0] = df_means_sem.loc['2020-05_COVID', 0]
        covId_means_sem.loc['2020-06', 0] = df_means_sem.loc['2020-06_COVID', 0]

        non_covId_means_avg = pd.DataFrame(index=df_means_avg_history.index)
        non_covId_means_sem = pd.DataFrame(index=df_means_sem_history.index)

        non_covId_means_avg.loc['2020-01', 0] = df_means_avg.loc['2020-01_Non_COVID', 0]
        non_covId_means_avg.loc['2020-02', 0] = df_means_avg.loc['2020-02_Non_COVID', 0]
        non_covId_means_avg.loc['2020-03', 0] = df_means_avg.loc['2020-03_Non_COVID', 0]
        non_covId_means_avg.loc['2020-04', 0] = df_means_avg.loc['2020-04_Non_COVID', 0]
        non_covId_means_avg.loc['2020-05', 0] = df_means_avg.loc['2020-05_Non_COVID', 0]
        non_covId_means_avg.loc['2020-06', 0] = df_means_avg.loc['2020-06_Non_COVID', 0]

        non_covId_means_sem.loc['2020-01', 0] = df_means_sem.loc['2020-01_Non_COVID', 0]
        non_covId_means_sem.loc['2020-02', 0] = df_means_sem.loc['2020-02_Non_COVID', 0]
        non_covId_means_sem.loc['2020-03', 0] = df_means_sem.loc['2020-03_Non_COVID', 0]
        non_covId_means_sem.loc['2020-04', 0] = df_means_sem.loc['2020-04_Non_COVID', 0]
        non_covId_means_sem.loc['2020-05', 0] = df_means_sem.loc['2020-05_Non_COVID', 0]
        non_covId_means_sem.loc['2020-06', 0] = df_means_sem.loc['2020-06_Non_COVID', 0]

        df_means_avg_history['COVID-19'] = covId_means_avg
        df_means_sem_history['COVID-19'] = covId_means_sem
        df_means_avg_history['Non COVID-19'] = non_covId_means_avg
        df_means_sem_history['Non COVID-19'] = non_covId_means_sem
        # all_counts_sum.loc[['2020-01_COVID','2020-02_COVID', '2020-03_COVID','2020-04_COVID','2020-05_COVID','2020-06_COVID'],:]
        vis.plt_acc_time_avgs_history(df_means_avg_history, df_means_sem_history,
                                      'COVID-19, Non COVID-19 and Total publications average time to acceptance')


    def get_acc_time_journals(self, df):
        metrics = utils.load_csv_data_to_df('journals_list_metrics_tmp.csv')
        metrics['SJR'] = pd.to_numeric(metrics['SJR'], errors='coerce')
        metrics['CiteScore'] = pd.to_numeric(metrics['CiteScore'], errors='coerce')
        metrics['SNIP'] = pd.to_numeric(metrics['SNIP'], errors='coerce')

        df_means = df[df['pub_name'].str.match('.*_mean')].copy()
        df_means['pub_name'] = df_means.apply(
            lambda row: pd.Series(row['pub_name'].replace('_acc_mean', '')), axis=1)

        df_means['SJR'] = df_means.apply(
            lambda row: pd.Series(self.get_SJR(row['pub_name'], metrics)), axis=1)

        df_means = df_means.sort_values('SJR', ascending=True, na_position='last')

        df_means.set_index('pub_name', inplace=True)
        df_means.rename(index={0: 'pub_name'}, inplace=True)
        df_means = df_means.drop(['New Scientist'])
        df_means = df_means.apply(pd.to_numeric, errors='coerce')

        df_means_hist = df_means.T.loc['2016-01':'2020-01']
        df_means_hist = df_means_hist.append(
            df_means.T.loc[['2020-02', '2020-03', '2020-04', '2020-05', '2020-06']])
        df_means_avg=df_means_hist.mean()
        cov_non_cov_means=df_means.T.loc[['2020-01_COVID','2020-01_Non_COVID',
                                               '2020-02_COVID','2020-02_Non_COVID',
                                               '2020-03_COVID','2020-03_Non_COVID',
                                               '2020-04_COVID','2020-04_Non_COVID',
                                               '2020-05_COVID','2020-05_Non_COVID',
                                               '2020-06_COVID','2020-06_Non_COVID']]


        print(cov_non_cov_means)
        title = 'COVID-19 vs Non COVID-19 Average time to acceptance in selected journals, compared with 2016-2020'
        vis.plt_covid_by_non_covid_acc_time(cov_non_cov_means.T, df_means_avg, title)


    def get_SJR(self,pub_name,metrics):
        pub_metrics = metrics[metrics['journal_name'] == pub_name]
        SJR = float(pub_metrics['SJR'])
        return(SJR)

    def get_top_publishing_countries(self,countries_dict):
        countries_dict_2020=countries_dict['2020']
        month_num=1
        month_dict=countries_dict_2020[str(month_num)]
        month_df=pd.DataFrame.from_dict(month_dict)
        month_df=month_df.T
        month_df = month_df.reset_index()
        month_df=month_df.sort_values(by='Country Total Articles',ascending=False, ignore_index=True)
        total_articles=month_df[0:1]['Country Total Articles'].values
        top_publishing_countries=month_df[1:6]
        top_publishing_countries.set_index('index', inplace=True)
        top_publishing_countries=top_publishing_countries['Country Total Articles']

        month_covid_dict = countries_dict_2020[str(month_num)+'_COVID']
        month_covid_df = pd.DataFrame.from_dict(month_covid_dict)
        month_covid_df = month_covid_df.T
        month_covid_df = month_covid_df.reset_index()
        month_covid_df = month_covid_df.sort_values(by='Country Total Articles', ascending=False)
        total_covid_articles = month_covid_df[0:1]['Country Total Articles'].values
        top_covid_publishing_countries = month_covid_df[1:6]
        top_covid_publishing_countries.set_index('index', inplace=True)
        top_covid_publishing_countries=top_covid_publishing_countries['Country Total Articles']


        print(month_df)


if __name__ == '__main__':
    start_time = datetime.now()
    print(start_time)
    utils=Utils(path=path)
    sdd=SciencedirectData()
    aaa=AuthorsAndCountries()
    at=AcceptanceTime()
    vis=Visualization()
    countries_dict=utils.load_obj("countries_dict")
    sdd.get_top_publishing_countries(countries_dict)
    # countries_dict=sdd.extract_countries_data()
    # utils.save_obj(countries_dict,"countries_dict")
    exit(0)
    df= sdd.get_journals_counts_by_year_and_month()
    utils.write_to_csv(df, "journals_counts_mean.csv")
    df=utils.load_csv_data_to_df("journals_counts_mean.csv")
    sdd.create_counts_table(df)
    sdd.get_acc_time_journals(df)
    exit(0)
    df = pd.DataFrame()
    name = 'Journal_of_Infection'
    name1 = 'Journal_of_Infection_history.csv'
    df=utils.load_csv_data_to_df(name1)
    name2=name='Journal_of_Infection_COVID.csv'
    df2 = utils.load_csv_data_to_df(name2)
    sdd.test_covid_in_data(df, df2)
    exit(0)
    name='Journal_of_Infection_COVID'
    obj_path = os.path.join(path, "ScienceDirect_citations_" + name + ".bib")
    papers = sdd.get_data(obj_path)
    utils.write_to_csv(papers, name + '.csv')
    utils.save_obj(papers, name)
    name = 'The_Lancet_Infectious_Diseases_history'
    obj_path = os.path.join(path, "ScienceDirect_citations_" + name + ".bib")
    papers = sdd.get_data(obj_path)
    utils.write_to_csv(papers, name + '.csv')
    utils.save_obj(papers, name)
    # for year in range(2016,2021):
    #     str_year=str(year)
    #     name='The_Lancet_Respiratory_Medicine_'
    #     obj_path = os.path.join(path, "ScienceDirect_citations_"+name+str_year+".bib")
    #     papers=pa.get_data(obj_path)
    #     df=df.append(papers)
    # utils.write_to_csv(papers,name+'history.csv')
    # utils.save_obj(papers,name+'history')
    end_time = datetime.now()
    print(end_time - start_time)
