from elsapy.elsclient import ElsClient
from elsapy.elsprofile import ElsAuthor, ElsAffil
from elsapy.elsdoc import FullDoc, AbsDoc
from elsapy.elssearch import ElsSearch
import json
from utils import Utils
from bio import Entrez
import pycountry
from parse_authors_special import ParseSpecialAuthors
import datetime
import pandas as pd
import itertools
from analyse_acceptance import AcceptanceAnalysis
from parse_scopus_csv import ParseScopusCsv
from sciencedirect_data import SciencedirectData
from parse_arxivs import ParseArxivs
import numpy as np
import imgkit
from visualization import Visualization

top_score_journals_names=['The Lancet', 'The Lancet Infectious Diseases', 'The Lancet Respiratory Medicine','The Lancet Global Health', 'The Lancet Public Health', 'Gastroenterology', 'Journal of the American Academy of Dermatology', 'Journal of Vascular Surgery' ]
selected_journals_names=['The Lancet Infectious Diseases', 'The Lancet Respiratory Medicine', 'The Lancet Global Health', 'Journal of Hospital Infection','International Journal of Infectious Diseases', 'Travel Medicine and Infectious Disease', 'European Urology', 'Psychiatry Research', 'Medical Hypotheses']
production_path=''
# local_path='D:\\shir\\study\\covid_19'
# local_path='D:\\shir\\study\\covid_19\\scopus'
local_path='D:\\shir\\study\\covid_19\\scopus\scienceDirectData'


def get_covid_growth(all_counts_sum):
    df = pd.DataFrame(
        columns=['Pub_name', 'SJR', 'COVID_Jan', 'Total_Jan', 'COVID_Feb', 'Total_Feb', 'COVID_Mar', 'Total_Mar',
                 'COVID_Apr', 'Total_Apr'])
    # obj1=utils.load_obj('countries_collab_'+str(month)+'_'+str(year))
    papers_2020 = utils.load_csv_data_to_df('2020_acceptance_time.csv')
    papers_covid = utils.load_csv_data_to_df('covid_acceptance_time.csv')
    papers_covid_jan = papers_covid[0:8]
    papers_covid_feb = papers_covid[8:38]
    papers_covid_mar = papers_covid[38:68]
    papers_covid_apr = papers_covid[68:]
    metrics = utils.load_csv_data_to_df('journals_list_metrics.csv')
    metrics['SJR'] = pd.to_numeric(metrics['SJR'], errors='coerce')
    metrics['CiteScore'] = pd.to_numeric(metrics['CiteScore'], errors='coerce')
    metrics['SNIP'] = pd.to_numeric(metrics['SNIP'], errors='coerce')

    metrics = metrics.sort_values(by='SJR')
    for row in papers_2020.iterrows():
        pub_name = row[1].Pub_name
        if pub_name == 'New Scientist' or pub_name=='Advances in Radiation Oncology' or pub_name=='Journal of Clinical Virology'\
                or pub_name=='Asian Journal of Psychiatry':
            continue
        # print(pub_name)
        pub_data_jan = papers_covid_jan['Pub_name'] == pub_name
        pub_data_feb = papers_covid_feb['Pub_name'] == pub_name
        pub_data_mar = papers_covid_mar['Pub_name'] == pub_name
        pub_data_apr = papers_covid_apr['Pub_name'] == pub_name
        pub_metrics = metrics[metrics['journal_name'] == pub_name]
        try:
            SJR = float(pub_metrics['SJR'])
        except:
            print('no SJR for {} '.format(pub_name))
            SJR = 0
        try:
            Total_Jan = row[1].count_nonzero_20201
            COVID_Jan = papers_covid_jan.loc[pub_data_jan, 'count_nonzero'].reset_index(drop=True)[0]
        except:
            COVID_Jan = 0

        try:
            Total_Feb = row[1].count_nonzero_20202
            COVID_Feb = papers_covid_feb.loc[pub_data_feb, 'count_nonzero'].reset_index(drop=True)[0]
        except:
            COVID_Feb = 0
        try:
            Total_Mar = row[1].count_nonzero_20203
            COVID_Mar = papers_covid_mar.loc[pub_data_mar, 'count_nonzero'].reset_index(drop=True)[0]
        except:
            COVID_Mar = 0
        try:
            Total_Apr = row[1].count_nonzero_20204
            COVID_Apr = papers_covid_apr.loc[pub_data_apr, 'count_nonzero'].reset_index(drop=True)[0]
        except:
            COVID_Apr = 0
        if COVID_Apr == 0:
            continue
        df = df.append(
            {'Pub_name': pub_name, 'SJR': SJR, 'COVID_Jan': COVID_Jan, 'Total_Jan': Total_Jan, 'COVID_Feb': COVID_Feb,
             'Total_Feb': Total_Feb, 'COVID_Mar': COVID_Mar, 'Total_Mar': Total_Mar, 'COVID_Apr': COVID_Apr,
             'Total_Apr': Total_Apr}, ignore_index=True)


    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    # get_longtitudal_avg_num_papers_change()

    df=df.set_index('Pub_name')
    df=df.apply(pd.to_numeric, errors='coerce')

    df_sums=df.sum()
    df_sums=df_sums.T
    df_sums=df_sums.drop(['SJR'])
    COV=dict()
    COV['January 2016'] = 0
    COV['February 2016'] = 0
    COV['March 2016'] = 0
    COV['April 2016'] = 0
    COV['January 2017'] = 0
    COV['February 2017'] = 0
    COV['March 2017'] = 0
    COV['April 2017'] = 0
    COV['January 2018'] = 0
    COV['February 2018'] = 0
    COV['March 2018'] = 0
    COV['April 2018'] = 0
    COV['January 2019'] = 0
    COV['February 2019'] = 0
    COV['March 2019'] = 0
    COV['April 2019'] = 0
    COV['January 2020'] = df_sums['COVID_Jan']
    COV['February 2020'] = df_sums['COVID_Feb']
    COV['March 2020'] = df_sums['COVID_Mar']
    COV['April 2020'] = df_sums['COVID_Apr']
    Totals = dict()
    Totals['January 2016'] = all_counts_sum['January 2016']
    Totals['February 2016'] = all_counts_sum['February 2016']
    Totals['March 2016'] = all_counts_sum['March 2016']
    Totals['April 2016'] = all_counts_sum['April 2016']
    Totals['January 2017'] = all_counts_sum['January 2017']
    Totals['February 2017'] = all_counts_sum['February 2017']
    Totals['March 2017'] = all_counts_sum['March 2017']
    Totals['April 2017'] = all_counts_sum['April 2017']
    Totals['January 2018'] = all_counts_sum['January 2018']
    Totals['February 2018'] = all_counts_sum['February 2018']
    Totals['March 2018'] = all_counts_sum['March 2018']
    Totals['April 2018'] = all_counts_sum['April 2018']
    Totals['January 2019'] = all_counts_sum['January 2019']
    Totals['February 2019'] = all_counts_sum['February 2019']
    Totals['March 2019'] = all_counts_sum['March 2019']
    Totals['April 2019'] = all_counts_sum['April 2019']
    Totals['January 2020']=df_sums['Total_Jan']
    Totals['February 2020']=df_sums['Total_Feb']
    Totals['March 2020']=df_sums['Total_Mar']
    Totals['April 2020']=df_sums['Total_Apr']
    vis.plt_journals_sums_old(COV,Totals, "COVID-19 and Total publications, Scopus journals")

    df['Percent_Jan'] = df.apply(lambda row: '{:.2%}'.format(100 * row.COVID_Jan / row.Total_Jan), axis=1)
    df['Percent_Feb'] = df.apply(lambda row: '{:.2%}'.format(100 * row.COVID_Feb / row.Total_Feb), axis=1)
    df['Percent_Mar'] = df.apply(lambda row: '{:.2%}'.format(100 * row.COVID_Mar / row.Total_Mar), axis=1)
    df['Percent_Apr'] = df.apply(lambda row: '{:.2%}'.format(100 * row.COVID_Apr / row.Total_Apr), axis=1)

    df=df.reindex(columns=['SJR', 'COVID_Jan', 'Total_Jan', 'Percent_Jan', 'COVID_Feb',
             'Total_Feb', 'Percent_Feb', 'COVID_Mar', 'Total_Mar', 'Percent_Mar', 'COVID_Apr',
             'Total_Apr', 'Percent_Apr'])

    df=df.sort_values('SJR', ascending=False, na_position='last')
    df_no_metrics=df.drop(['SJR'], axis=1)

    # vis.plt_covid_by_totals_journals(df_no_metrics)
    # df_no_metrics=df_no_metrics.T
    # vis.plt_journals_counts(df_no_metrics)
    # csv = df.to_csv()

    # df.style.set_table_styles([{'selector': '', 'props': [('border', '4px solid #7a7')]}])
    # tex=df.to_latex()
    # text_file = open("covid_growth_table.tex", "w")
    # text_file.write(tex)
    # text_file.close()

    print(df)

def get_non_covid_time_to_acceptance(num_covid_papers,covid_acc_time=1,num_total_papers=1,total_acc_time=1):
    total_total_time=num_total_papers*total_acc_time
    covid_total_time=num_covid_papers*covid_acc_time
    if covid_acc_time>0 and num_covid_papers==0:
        print('mismatch- acc time > 0 but num papers =0')
    if covid_total_time==0:
        num_covid_papers=0
    try:
        non_covid_acc_time=(total_total_time-covid_total_time)/(num_total_papers-num_covid_papers)
    except:
        non_covid_acc_time=total_acc_time
    return non_covid_acc_time


def fix_data_for_acc_time_avg(df):
    df.loc[df['COVID_acc_time_mean_Jan']==0, 'Total_acc_time_meanJan']=0
    df.loc[df['COVID_acc_time_mean_Jan'] == 0, 'Total_papers_jan']=0
    df.loc[df['COVID_acc_time_mean_Jan'] == 0, 'non_COVID_acc_time_mean_Jan'] = 0
    df.loc[df['COVID_acc_time_mean_Feb'] == 0, 'Total_acc_time_meanFeb'] = 0
    df.loc[df['COVID_acc_time_mean_Feb'] == 0, 'Total_papers_feb']=0
    df.loc[df['COVID_acc_time_mean_Feb'] == 0, 'non_COVID_acc_time_mean_Feb'] = 0
    df.loc[df['COVID_acc_time_mean_Mar'] == 0, 'Total_acc_time_meanMar'] = 0
    df.loc[df['COVID_acc_time_mean_Mar'] == 0, 'Total_papers_mar']=0
    df.loc[df['COVID_acc_time_mean_Mar'] == 0, 'non_COVID_acc_time_mean_Mar'] = 0
    df.loc[df['COVID_acc_time_mean_Apr'] == 0, 'Total_acc_time_meanApr'] = 0
    df.loc[df['COVID_acc_time_mean_Apr'] == 0, 'Total_papers_mar']=0
    df.loc[df['COVID_acc_time_mean_Apr'] == 0, 'non_COVID_acc_time_mean_Apr'] = 0


    return df


def get_COVID_avg_time_to_acceptance(all_acc_time_means_avg, avgs_for_journals_history):

    df = pd.DataFrame(
    columns=['Pub_name', 'COVID_acc_time_mean_Jan', 'Total_acc_time_meanJan', 'COVID_acc_time_mean_Feb', 'Total_acc_time_meanFeb', 'COVID_acc_time_mean_Mar', 'Total_acc_time_meanMar',
                     'COVID_acc_time_mean_Apr', 'Total_acc_time_meanApr'])
    # obj1=utils.load_obj('countries_collab_'+str(month)+'_'+str(year))
    papers_2020 = utils.load_csv_data_to_df('2020_acceptance_time.csv')
    papers_covid = utils.load_csv_data_to_df('covid_acceptance_time.csv')
    papers_covid_jan = papers_covid[0:8]
    papers_covid_feb = papers_covid[8:38]
    papers_covid_mar = papers_covid[38:68]
    papers_covid_apr = papers_covid[68:]
    for row in papers_2020.iterrows():
        pub_name = row[1].Pub_name
        if pub_name == 'New Scientist':
            continue
        # print(pub_name)
        pub_data_jan = papers_covid_jan['Pub_name'] == pub_name
        pub_data_feb = papers_covid_feb['Pub_name'] == pub_name
        pub_data_mar = papers_covid_mar['Pub_name'] == pub_name
        pub_data_apr = papers_covid_apr['Pub_name'] == pub_name

        non_COVID_acc_time_mean_Jan, non_COVID_acc_time_mean_Feb, non_COVID_acc_time_mean_Mar, non_COVID_acc_time_mean_Apr = '0','0','0','0'
        Total_acc_time_meanJan, Total_papers_jan, COVID_acc_time_mean_Jan, covid_papers_jan='0','0','0','0'
        Total_acc_time_meanFeb, Total_papers_feb, COVID_acc_time_mean_Feb, covid_papers_feb = '0', '0', '0', '0'
        Total_acc_time_meanMar, Total_papers_mar, COVID_acc_time_mean_Mar, covid_papers_mar = '0', '0', '0', '0'
        Total_acc_time_meanApr, Total_papers_apr, COVID_acc_time_mean_Apr, covid_papers_apr = '0', '0', '0', '0'
        try:
            Total_acc_time_meanJan = row[1].mean_20201
            Total_papers_jan=row[1].count_nonzero_20201
            # total_time_jan=Total_papers_jan*Total_acc_time_meanJan
            COVID_acc_time_mean_Jan = papers_covid_jan.loc[pub_data_jan, 'mean'].reset_index(drop=True)[0]
            covid_papers_jan=papers_covid_jan.loc[pub_data_jan, 'count_nonzero'].reset_index(drop=True)[0]
            # covid_total_time_jan=covid_papers_jan*COVID_acc_time_mean_Jan
            # non_COVID_acc_time_mean_Jan=(total_time_jan-covid_total_time_jan)/(Total_papers_jan-covid_papers_jan)
        except:
            COVID_acc_time_mean_Jan = 0
            non_COVID_acc_time_mean_Jan=Total_acc_time_meanJan

        try:
            Total_acc_time_meanFeb = row[1].mean_20202
            Total_papers_feb=row[1].count_nonzero_20202
            # total_time_feb=Total_papers_feb*Total_acc_time_meanFeb
            COVID_acc_time_mean_Feb = papers_covid_feb.loc[pub_data_feb, 'mean'].reset_index(drop=True)[0]
            covid_papers_feb = papers_covid_feb.loc[pub_data_feb, 'count_nonzero'].reset_index(drop=True)[0]
            # covid_total_time_feb = covid_papers_feb * COVID_acc_time_mean_Feb
            # non_COVID_acc_time_mean_Feb = (total_time_feb - covid_total_time_feb) / (
            #             Total_papers_feb - covid_papers_feb)
        except:
            COVID_acc_time_mean_Feb = 0
            non_COVID_acc_time_mean_Feb=Total_acc_time_meanFeb

        try:
            Total_acc_time_meanMar = row[1].mean_20203
            Total_papers_mar=row[1].count_nonzero_20203
            # total_time_mar=Total_papers_mar*Total_acc_time_meanMar
            COVID_acc_time_mean_Mar = papers_covid_mar.loc[pub_data_mar, 'mean'].reset_index(drop=True)[0]
            covid_papers_mar = papers_covid_mar.loc[pub_data_mar, 'count_nonzero'].reset_index(drop=True)[0]
            # covid_total_time_mar = covid_papers_mar * COVID_acc_time_mean_Mar
            # non_COVID_acc_time_mean_Mar = (total_time_mar - covid_total_time_mar) / (
            #         Total_papers_mar - covid_papers_mar)
        except:
            COVID_acc_time_mean_Mar = 0
            non_COVID_acc_time_mean_Mar=Total_acc_time_meanMar

        try:
            Total_acc_time_meanApr = row[1].mean_20204
            Total_papers_apr=row[1].count_nonzero_20204
            # total_time_apr=Total_papers_apr*Total_acc_time_meanApr
            COVID_acc_time_mean_Apr = papers_covid_apr.loc[pub_data_apr, 'mean'].reset_index(drop=True)[0]
            covid_papers_apr = papers_covid_apr.loc[pub_data_apr, 'count_nonzero'].reset_index(drop=True)[0]
            # covid_total_time_apr = covid_papers_apr * COVID_acc_time_mean_Apr
            # non_COVID_acc_time_mean_Apr = (total_time_apr - covid_total_time_apr) / (
            #         Total_papers_apr - covid_papers_apr)
        except:
            COVID_acc_time_mean_Apr = 0
            non_COVID_acc_time_mean_Apr=Total_acc_time_meanApr


        if COVID_acc_time_mean_Apr == 0:
            continue
        df = df.append(
            {'Pub_name': pub_name,  'COVID_acc_time_mean_Jan': COVID_acc_time_mean_Jan, 'covid_papers_jan': covid_papers_jan,
             'Total_acc_time_meanJan': Total_acc_time_meanJan, 'Total_papers_jan':Total_papers_jan, 'non_COVID_acc_time_mean_Jan': non_COVID_acc_time_mean_Jan,
             'COVID_acc_time_mean_Feb': COVID_acc_time_mean_Feb, 'covid_papers_feb':covid_papers_feb,
             'Total_acc_time_meanFeb': Total_acc_time_meanFeb, 'Total_papers_feb': Total_papers_feb, 'non_COVID_acc_time_mean_Feb': non_COVID_acc_time_mean_Feb,
             'COVID_acc_time_mean_Mar': COVID_acc_time_mean_Mar, 'covid_papers_mar':covid_papers_mar,
             'Total_acc_time_meanMar': Total_acc_time_meanMar, 'Total_papers_mar':Total_papers_mar, 'non_COVID_acc_time_mean_Mar': non_COVID_acc_time_mean_Mar,
             'COVID_acc_time_mean_Apr': COVID_acc_time_mean_Apr, 'covid_papers_apr': covid_papers_apr,
             'Total_acc_time_meanApr': Total_acc_time_meanApr, 'Total_papers_apr': Total_papers_apr,'non_COVID_acc_time_mean_Apr': non_COVID_acc_time_mean_Apr},
             ignore_index=True)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    df=df.set_index('Pub_name')
    df=df.apply(pd.to_numeric, errors='coerce')
    df2=df.copy()
    df2=fix_data_for_acc_time_avg(df2)
    df[['non_COVID_acc_time_mean_Jan']]=df.apply(
        lambda row: pd.Series(get_non_covid_time_to_acceptance(row['covid_papers_jan'], row['COVID_acc_time_mean_Jan'],
                                                               row['Total_papers_jan'], row['Total_acc_time_meanJan'])), axis=1)

    df[['non_COVID_acc_time_mean_Feb']] = df.apply(
        lambda row: pd.Series(get_non_covid_time_to_acceptance(row['covid_papers_feb'], row['COVID_acc_time_mean_Feb'],
                                                               row['Total_papers_feb'], row['Total_acc_time_meanFeb'])),
        axis=1)

    df[['non_COVID_acc_time_mean_Mar']] = df.apply(
        lambda row: pd.Series(get_non_covid_time_to_acceptance(row['covid_papers_mar'], row['COVID_acc_time_mean_Mar'],
                                                               row['Total_papers_mar'], row['Total_acc_time_meanMar'])),
        axis=1)

    df[['non_COVID_acc_time_mean_Apr']] = df.apply(
        lambda row: pd.Series(get_non_covid_time_to_acceptance(row['covid_papers_apr'], row['COVID_acc_time_mean_Apr'],
                                                               row['Total_papers_apr'], row['Total_acc_time_meanApr'])),
        axis=1)

    selected_df=(df.T[selected_journals_names]).T
    selected_df=selected_df.iloc[::-1]
    avgs_for_journals_history=avgs_for_journals_history.iloc[::-1]
    title='Average time to acceptance in selected journals, in first four months of 2020 compared with 2016-2019'
    # vis.plt_covid_by_non_covid_acc_time_old(selected_df, avgs_for_journals_history, title)
    print('done covid by non covid, now averages over months')
    df2=df2.replace(0,np.nan)
    df_avg=df2.mean()
    df_sem=df2.sem()
    all_acc_time_means_avg=all_acc_time_means_avg.T
    all_acc_time_means_avg['January 2020'][0]=df_avg.T['Total_acc_time_meanJan']
    all_acc_time_means_avg['January 2020'][1] = df_sem.T['Total_acc_time_meanJan']
    all_acc_time_means_avg['February 2020'][0] = df_avg.T['Total_acc_time_meanFeb']
    all_acc_time_means_avg['February 2020'][1] = df_sem.T['Total_acc_time_meanFeb']
    all_acc_time_means_avg['March 2020'][0] = df_avg.T['Total_acc_time_meanMar']
    all_acc_time_means_avg['March 2020'][1] = df_sem.T['Total_acc_time_meanMar']
    all_acc_time_means_avg['April 2020'][0] = df_avg.T['Total_acc_time_meanApr']
    all_acc_time_means_avg['April 2020'][1] = df_sem.T['Total_acc_time_meanApr']

    all_acc_time_means_sem=pd.DataFrame()
    all_acc_time_means_sem[0]=all_acc_time_means_avg.loc[1].copy()
    all_acc_time_means_sem=all_acc_time_means_sem.T

    all_acc_time_means_avg.loc[1]=[np.nan,np.nan,np.nan,np.nan,
                       np.nan,np.nan,np.nan,np.nan,
                       np.nan,np.nan,np.nan,np.nan,
                       np.nan,np.nan,np.nan,np.nan,
                       df_avg.T['COVID_acc_time_mean_Jan'],df_avg.T['COVID_acc_time_mean_Feb'],df_avg.T['COVID_acc_time_mean_Mar'],df_avg.T['COVID_acc_time_mean_Apr']]
    all_acc_time_means_sem.loc[1] = [np.nan, np.nan, np.nan, np.nan,
                                     np.nan, np.nan, np.nan, np.nan,
                                     np.nan, np.nan, np.nan, np.nan,
                                     np.nan, np.nan, np.nan, np.nan,
                                     df_sem.T['COVID_acc_time_mean_Jan'], df_sem.T['COVID_acc_time_mean_Feb'],
                                     df_sem.T['COVID_acc_time_mean_Mar'], df_sem.T['COVID_acc_time_mean_Apr']]

    all_acc_time_means_avg.loc[2] = [np.nan, np.nan, np.nan, np.nan,
                         np.nan, np.nan, np.nan, np.nan,
                         np.nan, np.nan, np.nan, np.nan,
                         np.nan, np.nan, np.nan, np.nan,
                         df_avg.T['non_COVID_acc_time_mean_Jan'], df_avg.T['non_COVID_acc_time_mean_Feb'],
                         df_avg.T['non_COVID_acc_time_mean_Mar'], df_avg.T['non_COVID_acc_time_mean_Apr']]

    all_acc_time_means_sem.loc[2] = [np.nan, np.nan, np.nan, np.nan,
                                     np.nan, np.nan, np.nan, np.nan,
                                     np.nan, np.nan, np.nan, np.nan,
                                     np.nan, np.nan, np.nan, np.nan,
                                     df_sem.T['non_COVID_acc_time_mean_Jan'], df_sem.T['non_COVID_acc_time_mean_Feb'],
                                     df_sem.T['non_COVID_acc_time_mean_Mar'], df_sem.T['non_COVID_acc_time_mean_Apr']]

    all_acc_time_means_avg.rename(index={0: 'Total', 1: 'COVID-19', 2: 'Non COVID-19'}, inplace=True)
    all_acc_time_means_sem.rename(index={0: 'Total', 1: 'COVID-19', 2: 'Non COVID-19'}, inplace=True)

    all_acc_time_means_avg=all_acc_time_means_avg.T
    all_acc_time_means_sem=all_acc_time_means_sem.T
    vis.plt_acc_time_avgs_history(all_acc_time_means_avg, all_acc_time_means_sem,'COVID-19, Non COVID-19 and Total publications average time to acceptance')
    print(all_acc_time_means_avg)



def get_longtitudal_avg_time_to_acceptance_changes():
    papers_2020 = utils.load_csv_data_to_df('2020_acceptance_time.csv')
    papers_2019 = utils.load_csv_data_to_df('2019_acceptance_time.csv')
    papers_2018 = utils.load_csv_data_to_df('2018_acceptance_time.csv')
    papers_2017 = utils.load_csv_data_to_df('2017_acceptance_time.csv')
    papers_2016 = utils.load_csv_data_to_df('2016_acceptance_time.csv')
    papers_2020_acc_time_mean=papers_2020[['Pub_name','mean_20201','mean_20202','mean_20203','mean_20204']]
    papers_2020_acc_time_mean.rename(columns={'mean_20201': 'January 2020', 'mean_20202': 'February 2020', 'mean_20203': 'March 2020','mean_20204': 'April 2020'}, inplace=True)

    papers_2019_acc_time_mean=papers_2019[['Pub_name','mean_20191','mean_20192','mean_20193','mean_20194']]
    papers_2019_acc_time_mean.rename(columns={'mean_20191': 'January 2019', 'mean_20192': 'February 2019',
                                      'mean_20193': 'March 2019', 'mean_20194': 'April 2019'},
                             inplace=True)
    papers_2018_acc_time_mean=papers_2018[['Pub_name','mean_20181','mean_20182','mean_20183','mean_20184']]
    papers_2018_acc_time_mean.rename(columns={'mean_20181': 'January 2018', 'mean_20182': 'February 2018',
                                      'mean_20183': 'March 2018', 'mean_20184': 'April 2018'},
                             inplace=True)
    papers_2017_acc_time_mean=papers_2017[['Pub_name','mean_20171','mean_20172','mean_20173','mean_20174']]
    papers_2017_acc_time_mean.rename(columns={'mean_20171': 'January 2017', 'mean_20172': 'February 2017',
                                      'mean_20173': 'March 2017', 'mean_20174': 'April 2017'},
                             inplace=True)
    papers_2016_acc_time_mean=papers_2016[['Pub_name','mean_20161','mean_20162','mean_20163','mean_20164']]
    papers_2016_acc_time_mean.rename(columns={'mean_20161': 'January 2016', 'mean_20162': 'February 2016',
                                      'mean_20163': 'March 2016', 'mean_20164': 'April 2016'},
                             inplace=True)
    all_acc_time_means=pd.merge(pd.merge(pd.merge(pd.merge(papers_2016_acc_time_mean,papers_2017_acc_time_mean),papers_2018_acc_time_mean),papers_2019_acc_time_mean),papers_2020_acc_time_mean)
    all_acc_time_means=all_acc_time_means.set_index('Pub_name')
    all_acc_time_means=all_acc_time_means.apply(pd.to_numeric, errors='coerce')
    all_acc_time_means=all_acc_time_means.T
    all_acc_time_means=all_acc_time_means.drop(['New Scientist','Infection Genetics and Evolution','Journal of Microbiology Immunology and Infection', 'The Lancet Gastroenterology & Hepatology', 'Infectious Disease Modelling', 'Asian Journal of Psychiatry', 'Brain  Behavior  and Immunity','Diabetes & Metabolic Syndrome: Clinical Research & Reviews', 'Advances in Radiation Oncology','British Journal of Anaesthesia'], axis=1)
    all_acc_time_means.reset_index(level=0, inplace=True)
    all_acc_time_means.rename(columns={'index': 'months'}, inplace=True)


    all_acc_time_means_avg = pd.DataFrame()
    all_acc_time_means_avg[0]=all_acc_time_means.mean(axis=1)
    all_acc_time_means_avg[1] = all_acc_time_means.sem(axis=1)
    all_acc_time_means_avg.set_index(all_acc_time_means['months'], inplace=True)
    avgs_for_journals=all_acc_time_means[0:16].mean()
    avgs_for_journals=avgs_for_journals[selected_journals_names]

    get_COVID_avg_time_to_acceptance(all_acc_time_means_avg, avgs_for_journals)
    # all_acc_time_means_avg.set_index(level=0, inplace=True, drop=True)


    title = 'Mean time to acceptance in years 2016-2020'
    # vis.plt_journals_counts(all_acc_time_means_avg, title)

    month_selected_journals_names = selected_journals_names
    month_selected_journals_names.insert(0, 'months')
    all_acc_time_means_selected = all_acc_time_means[month_selected_journals_names]

    title='Mean time to acceptance in years 2016-2020'
    # vis.plt_journals_counts(all_acc_time_means_selected, title)

def get_longtitudal_num_papers_changes():
    papers_2020 = utils.load_csv_data_to_df('2020_acceptance_time.csv')
    papers_2019 = utils.load_csv_data_to_df('2019_acceptance_time.csv')
    papers_2018 = utils.load_csv_data_to_df('2018_acceptance_time.csv')
    papers_2017 = utils.load_csv_data_to_df('2017_acceptance_time.csv')
    papers_2016 = utils.load_csv_data_to_df('2016_acceptance_time.csv')
    papers_2020_count=papers_2020[['Pub_name','count_nonzero_20201','count_nonzero_20202','count_nonzero_20203','count_nonzero_20204']]
    papers_2020_count.rename(columns={'count_nonzero_20201': 'January 2020', 'count_nonzero_20202': 'February 2020', 'count_nonzero_20203': 'March 2020','count_nonzero_20204': 'April 2020'}, inplace=True)

    papers_2019_count=papers_2019[['Pub_name','count_nonzero_20191','count_nonzero_20192','count_nonzero_20193','count_nonzero_20194']]
    papers_2019_count.rename(columns={'count_nonzero_20191': 'January 2019', 'count_nonzero_20192': 'February 2019',
                                      'count_nonzero_20193': 'March 2019', 'count_nonzero_20194': 'April 2019'},
                             inplace=True)
    papers_2018_count=papers_2018[['Pub_name','count_nonzero_20181','count_nonzero_20182','count_nonzero_20183','count_nonzero_20184']]
    papers_2018_count.rename(columns={'count_nonzero_20181': 'January 2018', 'count_nonzero_20182': 'February 2018',
                                      'count_nonzero_20183': 'March 2018', 'count_nonzero_20184': 'April 2018'},
                             inplace=True)
    papers_2017_count=papers_2017[['Pub_name','count_nonzero_20171','count_nonzero_20172','count_nonzero_20173','count_nonzero_20174']]
    papers_2017_count.rename(columns={'count_nonzero_20171': 'January 2017', 'count_nonzero_20172': 'February 2017',
                                      'count_nonzero_20173': 'March 2017', 'count_nonzero_20174': 'April 2017'},
                             inplace=True)
    papers_2016_count=papers_2016[['Pub_name','count_nonzero_20161','count_nonzero_20162','count_nonzero_20163','count_nonzero_20164']]
    papers_2016_count.rename(columns={'count_nonzero_20161': 'January 2016', 'count_nonzero_20162': 'February 2016',
                                      'count_nonzero_20163': 'March 2016', 'count_nonzero_20164': 'April 2016'},
                             inplace=True)
    all_counts=pd.merge(pd.merge(pd.merge(pd.merge(papers_2016_count,papers_2017_count),papers_2018_count),papers_2019_count),papers_2020_count)
    all_counts=all_counts.set_index('Pub_name')
    all_counts=all_counts.apply(pd.to_numeric, errors='coerce')
    all_counts = all_counts.T
    all_counts=all_counts.drop(['Advances in Radiation Oncology','Journal of Clinical Virology'], axis=1)
    all_counts=all_counts.drop(['New Scientist','Infection Genetics and Evolution',
                                'Journal of Microbiology Immunology and Infection',
                                'The Lancet Gastroenterology & Hepatology',
                                'Infectious Disease Modelling', 'Asian Journal of Psychiatry', 'Brain  Behavior  and Immunity','Diabetes & Metabolic Syndrome: Clinical Research & Reviews'], axis=1)
    all_counts = all_counts.drop(
        ['Journal of Genetics and Genomics', 'The Lancet Psychiatry', 'Microbes and Infection',
         'International Journal of Antimicrobial Agents'], axis=1)
    all_counts=all_counts.T
    all_counts_sum = all_counts.sum()
    get_covid_growth(all_counts_sum)

    all_counts.reset_index(level=0, inplace=True)
    all_counts.rename(columns={'index': 'months'}, inplace=True)
    month_top_score_journals_names=top_score_journals_names
    month_top_score_journals_names.insert(0,'months')
    top_score_journals=all_counts[month_top_score_journals_names]
    title='Publication counts in years 2016-2020 in journals with highest SJR scores'
    vis.plt_journals_counts(top_score_journals, title)


def get_journals_publication_growth_plot():
    df = utils.load_csv_data_to_df("journals_counts_mean.csv")
    scd.get_history_counts(df,vis)

def get_arxivs_publication_growth_plot():
    data = utils.load_csv_data_to_df('covid_counts_preprints.csv')
    pa.get_longtitudal_num_papers_changes(data,vis)

def get_avg_time_to_acc_journals():
    df = utils.load_csv_data_to_df("journals_counts_mean.csv")
    scd.get_acc_time_journals(df, utils,vis)

def get_avg_time_to_acc_journals_history():
    df = utils.load_csv_data_to_df("journals_counts_mean.csv")
    scd.get_acc_time(df, vis)

def get_top_publishing_countries_history():
    countries_dict = utils.load_obj("countries_dict")
    scd.get_top_publishing_countries_history(countries_dict,vis)

def get_top_collab_countries_diversity():
    countries_dict = utils.load_obj("countries_dict")
    scd.get_countries_collab_covid(countries_dict,vis)

def get_top_collab_countries_diversity_history():
    countries_dict = utils.load_obj("countries_dict")
    scd.get_country_collab_diversity_history(countries_dict,vis)

def get_top_collab_countries_num_papers():
    countries_dict = utils.load_obj("countries_dict")
    scd.get_collab_by_num_papers_covid(countries_dict,vis)

def get_top_collab_countries_num_papers_history():
    countries_dict = utils.load_obj("countries_dict")
    scd.get_collab_by_num_papers_history(countries_dict,vis)


if __name__ == '__main__':

    # imgkitoptions = {"format": format}
    #
    # imgkit.from_file("covid_growth_table.html", "covid_growth_table", options=imgkitoptions)
    #
    # exit()
    print(datetime.datetime.now())
    utils=Utils(path=local_path)
    aa=AcceptanceAnalysis()
    psc=ParseScopusCsv()
    scd = SciencedirectData()
    pa=ParseArxivs()
    vis=Visualization()
    # get_arxivs_publication_growth_plot()
    # get_journals_publication_growth_plot()
    get_avg_time_to_acc_journals()
    # get_avg_time_to_acc_journals_history()
    # get_top_collab_countries_diversity()
    # get_top_collab_countries_diversity_history()
    # get_top_collab_countries_num_papers()
    get_top_collab_countries_num_papers_history()
    exit(0)
    # month = 1
    # year=2020
    # get_longtitudal_num_papers_changes()
    get_longtitudal_avg_time_to_acceptance_changes()
    # get_COVID_avg_time_to_acceptance()
    # get_covid_growth()
    # arrays=[np.hstack([['January']*2,['February']*2,['March']*2,['April']*2]),['COVID','Total']*4]
    # columns=pd.MultiIndex.from_arrays(arrays, names=['Months',['Journal']])

    # for row in papers.iterrows():
    #     DOI=row[1].DOI
    #     res=aa.get_data_from_doi(DOI)
    # papers[['Date_Received', 'Date_Accepted', 'Date_online', 'Acceptance_Time']] = papers.apply(
    #     lambda row: pd.Series(scd.dates_calc(aa.get_data_from_doi(row['DOI']))), axis=1)
    # a = papers.loc[papers['Date_online'].str.match('.*/08/2020'), :]
    # b = papers.loc[papers['Date_online'].str.match('.*/05/2020'), :]
    # c = papers.loc[papers['Date_online'].str.match('.*/06/2020'), :]
    # d = papers.loc[papers['Date_online'].str.match('.*/07/2020'), :]


    print(datetime.datetime.now())
