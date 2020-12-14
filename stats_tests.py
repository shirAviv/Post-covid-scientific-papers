
import matplotlib.pyplot as plt
import random
from scipy import stats
from sciencedirect_data import SciencedirectData
from utils import Utils
import numpy as np
import pandas as pd
from scipy.stats import shapiro
import math
path='D:\\shir\\study\\covid_19\\scopus\scienceDirectData'

class StatsTests:

    def temp_test(self):
        random.seed(20) #for results to be recreated

        N = 50 #number of samples to take from each population
        a = [random.gauss(55,20) for x in range(N)] #take N samples from population A
        b = [random.gauss(50,15) for x in range(N+10)] #take N samples from population B

        plt.plot(sorted(a))
        plt.plot(b)
        plt.title("Independent Sample T-Test")
        plt.show()

        tStat, pValue = stats.ttest_ind(a, b, equal_var = False) #run independent sample T-Test
        print("P-Value:{0} T-Statistic:{1}".format(pValue,tStat)) #print the P-Value and the T-Statistic
        ##P-Value:0.017485741540118758 T-Statistic:2.421942924642376

    def test_for_normal_dist(self, df):
        df=df.T
        df[['normal dist']] = df.apply(
            lambda row: pd.Series(self.is_normal_dist(row)), axis=1)
        # df['normal dist']=df.apply(lambda row)


    def is_normal_dist(vals):
            # data = [1, 1.2, 0.2, 0.3, -1, -0.2, -0.6, -0.8, 0.8, 0.1]
        stat, p = shapiro(vals)
        print('stat=%.3f, p=%.3f' % (stat, p))
        if p > 0.05:
            print("Data follows Normal Distribution")
            return True
        else:
            print("Data does not follow Normal Distribution")
            return False

    def run_ttest_covid_two_groups(self, acc_df, group_1, group_2, timedelta=True):
        if timedelta:
            a = acc_df.loc[group_1].apply(lambda row: row.days)
            b = acc_df.loc[group_2].apply(lambda row: row.days)
        else:
            a = acc_df.loc[group_1]
            b = acc_df.loc[group_2]
        a_mean = np.mean(np.ma.masked_invalid(a.values))
        a_std = np.std(np.ma.masked_invalid(a.values))
        b_mean = np.mean(np.ma.masked_invalid(b.values))
        b_std = np.std(np.ma.masked_invalid(b.values))

        tStat, pValue = stats.ttest_ind(b.values, a.values, equal_var=False, nan_policy='omit')
        pValue=pValue/2
        return a_mean, a_std, b_mean, b_std, pValue, tStat

    def run_wtest_covid_non_covid(self, acc_df, group_1, group_2, timedelta=True):
        if timedelta:
            a = acc_df.loc[group_1].apply(lambda row: row.days)
            b = acc_df.loc[group_2].apply(lambda row: row.days)
        else:
            a = acc_df.loc[group_1]
            b = acc_df.loc[group_2]
        a_mean = np.mean(np.ma.masked_invalid(a.values))
        a_std = np.std(np.ma.masked_invalid(a.values))
        b_mean = np.mean(np.ma.masked_invalid(b.values))
        b_std = np.std(np.ma.masked_invalid(b.values))

        wStat, pValue=stats.wilcoxon(b.values, a.values)
        # tStat, pValue = stats.ttest_ind(b.values, a.values, equal_var=False, nan_policy='omit')
        # pValue=pValue/2
        return a_mean, a_std, b_mean, b_std, pValue, wStat

    def run_ttest_with_groups(self, df, gr1_1, gr1_2, gr1_3, gr1_4, gr1_5, gr1_6, gr2_1, gr2_2, gr2_3, gr2_4, gr2_5, gr2_6, timedelta=True):
        if timedelta:
            a = df.loc[gr1_1].apply(lambda row: row.days)
            a = a.append(df.loc[gr1_2].apply(lambda row: row.days))
            a = a.append(df.loc[gr1_3].apply(lambda row: row.days))
            a = a.append(df.loc[gr1_4].apply(lambda row: row.days))
            a = a.append(df.loc[gr1_5].apply(lambda row: row.days))
            a = a.append(df.loc[gr1_6].apply(lambda row: row.days))

            b = df.loc[gr2_1].apply(lambda row: row.days)
            b = b.append(df.loc[gr2_2].apply(lambda row: row.days))
            b = b.append(df.loc[gr2_3].apply(lambda row: row.days))
            b = b.append(df.loc[gr2_4].apply(lambda row: row.days))
            b = b.append(df.loc[gr2_5].apply(lambda row: row.days))
            b = b.append(df.loc[gr2_6].apply(lambda row: row.days))
        else:
            df_with_zeros=df.fillna(0)
            a = df_with_zeros.loc[gr1_1] + df_with_zeros.loc[gr1_2] + df_with_zeros.loc[gr1_3] + df_with_zeros.loc[gr1_4] + df_with_zeros.loc[gr1_5] + df_with_zeros.loc[gr1_6]

            b = df_with_zeros.loc[gr2_1] + df_with_zeros.loc[gr2_2] + df_with_zeros.loc[gr2_3] + df_with_zeros.loc[gr2_4] + df_with_zeros.loc[gr2_5] + df_with_zeros.loc[gr2_6]
        a=a.replace(0,np.nan)
        b = b.replace(0, np.nan)
        a_mean = np.mean(np.ma.masked_invalid(a.values))
        a_std = np.std(np.ma.masked_invalid(a.values))
        b_mean = np.mean(np.ma.masked_invalid(b.values))
        b_std = np.std(np.ma.masked_invalid(b.values))

        tStat, pValue = stats.ttest_ind(b.values, a.values, equal_var=False, nan_policy='omit')
        pValue=pValue/2
        return a_mean, a_std, b_mean, b_std, pValue, tStat

    def run_wtest(self, acc_df, gr1_1, gr1_2, gr1_3, gr1_4, gr1_5, gr1_6, gr2_1, gr2_2, gr2_3, gr2_4, gr2_5, gr2_6, timedelta=True):
        if timedelta:
            a = acc_df.loc[gr1_1].apply(lambda row: row.days)
            a = a.append(acc_df.loc[gr1_2].apply(lambda row: row.days))
            a = a.append(acc_df.loc[gr1_3].apply(lambda row: row.days))
            a = a.append(acc_df.loc[gr1_4].apply(lambda row: row.days))
            a = a.append(acc_df.loc[gr1_5].apply(lambda row: row.days))
            a = a.append(acc_df.loc[gr1_6].apply(lambda row: row.days))

            b = acc_df.loc[gr2_1].apply(lambda row: row.days)
            b = b.append(acc_df.loc[gr2_2].apply(lambda row: row.days))
            b = b.append(acc_df.loc[gr2_3].apply(lambda row: row.days))
            b = b.append(acc_df.loc[gr2_4].apply(lambda row: row.days))
            b = b.append(acc_df.loc[gr2_5].apply(lambda row: row.days))
            b = b.append(acc_df.loc[gr2_6].apply(lambda row: row.days))
        else:
            a = acc_df.loc[gr1_1]
            a = a.append(acc_df.loc[gr1_2])
            a = a.append(acc_df.loc[gr1_3])
            a = a.append(acc_df.loc[gr1_4])
            a = a.append(acc_df.loc[gr1_5])
            a = a.append(acc_df.loc[gr1_6])

            b = acc_df.loc[gr2_1]
            b = b.append(acc_df.loc[gr2_2])
            b = b.append(acc_df.loc[gr2_3])
            b = b.append(acc_df.loc[gr2_4])
            b = b.append(acc_df.loc[gr2_5])
            b = b.append(acc_df.loc[gr2_6])

        a_mean = np.mean(np.ma.masked_invalid(a.values))
        a_std = np.std(np.ma.masked_invalid(a.values))
        b_mean = np.mean(np.ma.masked_invalid(b.values))
        b_std = np.std(np.ma.masked_invalid(b.values))

        wStat, pValue = stats.wilcoxon(b.values, a.values)
        pValue=pValue
        return a_mean, a_std, b_mean, b_std, pValue, wStat

    def run_acc_ttests(self):
        acc_df = utils.load_obj('acceptance_df')
        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(acc_df, '2020_02_non_COVID', '2020_02_COVID')
        print("Feb -Non COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Feb -COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Feb -P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(acc_df, '2020_03_non_COVID', '2020_03_COVID')
        print("Mar -Non COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Mar -COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Mar -P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(acc_df, '2020_04_non_COVID', '2020_04_COVID')
        print("Apr -Non COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Apr -COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Apr- P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(acc_df, '2020_05_non_COVID', '2020_05_COVID')
        print("May -Non COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("May -COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("May -P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(acc_df, '2020_06_non_COVID', '2020_06_COVID')
        print("June -Non COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("June -COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("June- P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(acc_df,
                                                                   '2020_01_non_COVID', '2020_02_non_COVID',
                                                                   '2020_03_non_COVID', '2020_04_non_COVID',
                                                                   '2020_05_non_COVID', '2020_06_non_COVID',
                                                                   '2020_01_COVID', '2020_02_COVID',
                                                                   '2020_03_COVID', '2020_04_COVID',
                                                                   '2020_05_COVID', '2020_06_COVID'
                                                                               )

        print("Total 2020 -Non COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2020 -COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2020- P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        # history
        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(acc_df, '2016_01', '2016_02',
                                                                   '2016_03', '2016_04',
                                                                   '2016_05', '2016_06',
                                                                   '2017_01', '2017_02',
                                                                   '2017_03', '2017_04',
                                                                   '2017_05', '2017_06'
                                                                               )

        print("Total 2016 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2017 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2016, 2017- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(acc_df,
            '2017_01', '2017_02',
            '2017_03', '2017_04',
            '2017_05', '2017_06',
            '2018_01', '2018_02',
            '2018_03', '2018_04',
            '2018_05', '2018_06',
                                                                               )

        print("Total 2017 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2018 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2017, 2018- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(acc_df,
            '2018_01', '2018_02',
            '2018_03', '2018_04',
            '2018_05', '2018_06',
            '2019_01', '2019_02',
            '2019_03', '2019_04',
            '2019_05', '2019_06',
                                                                               )

        print("Total 2018 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2019 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2018, 2019- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(acc_df,
            '2019_01', '2019_02',
            '2019_03', '2019_04',
            '2019_05', '2019_06',
            '2020_01', '2020_02',
            '2020_03', '2020_04',
            '2020_05', '2020_06',
                                                                               )

        print("Total 2019 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2020 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2019, 2020- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(acc_df,
                                                                   '2019_01', '2019_02',
                                                                   '2019_03', '2019_04',
                                                                   '2019_05', '2019_06',
                                                                   '2020_01_non_COVID', '2020_02_non_COVID',
                                                                   '2020_03_non_COVID', '2020_04_non_COVID',
                                                                   '2020_05_non_COVID', '2020_06_non_COVID',
                                                                               )

        print("Total 2019 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2020 non COVID - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2019, 2020 non COVID- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

    def run_collab_countries_ttests(self, df):

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df, '2020_02_COVID',
                                                                                   '2020_02_non_COVID', timedelta=False)
        print("Feb -COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Feb -Non COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Feb -P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df, '2020_03_COVID',
                                                                                   '2020_03_non_COVID', timedelta=False)
        print("Mar -COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Mar -Non COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Mar -P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df, '2020_04_COVID',
                                                                                   '2020_04_non_COVID', timedelta=False)
        print("Apr -COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Apr -Non COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Apr- P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df, '2020_05_COVID',
                                                                                   '2020_05_non_COVID', timedelta=False)
        print("May -COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("May -Non COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("May -P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df, '2020_06_COVID',
                                                                                   '2020_06_non_COVID', timedelta=False)
        print("June -COVID mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("June -Non COVID mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("June- P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df,
                                                                   '2020_non_COVID_count', '2020_COVID_count',
                                                                               timedelta=False)

        print("Total 2020 -Non COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2020 -COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2020- P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        # history
        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df,
                                                                   '2017_count', '2016_count',
                                                                               timedelta=False)

        print("Total 2017 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2016 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2016, 2017- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df,
                                                                   '2018_count', '2017_count',
                                                                               timedelta=False)

        print("Total 2018 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2017 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2017, 2018- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df,
                                                                   '2019_count', '2018_count',
                                                                               timedelta=False)

        print("Total 2019 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2018 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2018, 2019- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df,
                                                                   '2020_count', '2019_count',
                                                                               timedelta=False)

        print("Total 2020 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2019 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2019, 2020- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df,
                                                                                     '2019_count','2020_COVID_count',
                                                                                    timedelta=False)

        print("Total 2019 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2020 COVID - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2019, 2020 COVID- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df,
                                                                                    '2020_non_COVID_count', '2019_count',
                                                                                    timedelta=False)

        print("Total 2020 non COVID - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2019 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2019, 2020 non COVID- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        # print(df)

    def run_collab_countries_by_num_papers_ttests(self, df):

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df, '2020_02_COVID',
                                                                                   '2020_02_non_COVID', timedelta=False)
        print("Feb -COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Feb -Non COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Feb -P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df, '2020_03_COVID',
                                                                                   '2020_03_non_COVID', timedelta=False)
        print("Mar -COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Mar -Non COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Mar -P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df, '2020_04_COVID',
                                                                                   '2020_04_non_COVID', timedelta=False)
        print("Apr -COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Apr -Non COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Apr- P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df, '2020_05_COVID',
                                                                                   '2020_05_non_COVID', timedelta=False)
        print("May -COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("May -Non COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("May -P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_covid_two_groups(df, '2020_06_COVID',
                                                                                   '2020_06_non_COVID', timedelta=False)
        print("June -COVID mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("June -Non COVID mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("June- P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(df,'2020_01_non_COVID', '2020_02_non_COVID',
                                                                   '2020_03_non_COVID', '2020_04_non_COVID',
                                                                   '2020_05_non_COVID', '2020_06_non_COVID',
                                                                   '2020_01_COVID', '2020_02_COVID',
                                                                   '2020_03_COVID', '2020_04_COVID',
                                                                   '2020_05_COVID', '2020_06_COVID',
                                                                               timedelta=False)

        print("Total 2020 -Non COVID-19 mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2020 -COVID-19 mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2020- P-Value:{0} T-Statistic:{1}".format(pValue, tStat))  # print the P-Value and the T-Statistic

        # history
        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(df,'2017_01', '2017_02',
                                                                   '2017_03', '2017_04',
                                                                   '2017_05', '2017_06',
                                                                   '2016_01', '2016_02',
                                                                   '2016_03', '2016_04',
                                                                   '2016_05', '2016_06',
                                                                               timedelta=False)

        print("Total 2017 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2016 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2016, 2017- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(df,'2018_01', '2018_02',
                                                                   '2018_03', '2018_04',
                                                                   '2018_05', '2018_06',
                                                                   '2017_01', '2017_02',
                                                                   '2017_03', '2017_04',
                                                                   '2017_05', '2017_06',
                                                                               timedelta=False)

        print("Total 2018 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2017 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2017, 2018- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(df,'2019_01', '2019_02',
                                                                   '2019_03', '2019_04',
                                                                   '2019_05', '2019_06',
                                                                   '2018_01', '2018_02',
                                                                   '2018_03', '2018_04',
                                                                   '2018_05', '2018_06',
                                                                               timedelta=False)

        print("Total 2019 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2018 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2018, 2019- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(df,
                                                                   '2020_01', '2020_02',
                                                                   '2020_03', '2020_04',
                                                                   '2020_05', '2020_06' ,
                                                                   '2019_01', '2019_02',
                                                                   '2019_03', '2019_04',
                                                                   '2019_05', '2019_06',
                                                                               timedelta=False)

        print("Total 2020 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2019 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2019, 2020- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(df,
                                                                   '2019_01', '2019_02',
                                                                   '2019_03', '2019_04',
                                                                   '2019_05', '2019_06',
                                                                   '2020_01_COVID', '2020_02_COVID',
                                                                   '2020_03_COVID', '2020_04_COVID',
                                                                   '2020_05_COVID', '2020_06_COVID',
                                                                                    timedelta=False)

        print("Total 2019 - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2020 COVID - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2019, 2020 COVID- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        a_mean, a_std, b_mean, b_std, pValue, tStat = st.run_ttest_with_groups(df,
                                                                   '2020_01_non_COVID', '2020_02_non_COVID',
                                                                   '2020_03_non_COVID', '2020_04_non_COVID',
                                                                   '2020_05_non_COVID', '2020_06_non_COVID',
                                                                   '2019_01', '2019_02',
                                                                   '2019_03', '2019_04',
                                                                   '2019_05', '2019_06',
                                                                                    timedelta=False)

        print("Total 2020 non COVID - mean:{0} std:{1}".format(a_mean, a_std))  # print the mean and std
        print("Total 2019 - mean:{0} std:{1}".format(b_mean, b_std))  # print the mean and std
        print("Total 2019, 2020 non COVID- P-Value:{0} T-Statistic:{1}".format(pValue,
                                                                     tStat))  # print the P-Value and the T-Statistic

        # print(df)

    def run_chi_single_vs_multi_country(self, df):
        current_df = df.loc[:, "2020_02_COVID":"2020_02_non_COVID"]
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("Feb -COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2, p))

        current_df = df.loc[:, "2020_03_COVID":"2020_03_non_COVID"]
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("Mar -COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2, p))

        current_df = df.loc[:, "2020_04_COVID":"2020_04_non_COVID"]
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("Apr -COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2, p))

        current_df = df.loc[:, "2020_05_COVID":"2020_05_non_COVID"]
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("May -COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2, p))

        current_df=df.loc[:,"2020_06_COVID" :"2020_06_non_COVID"]
        current_df=current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("June -COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2,p))

        current_df = df.loc[:, "2020_01_COVID":"2020_01_non_COVID"].copy()
        current_df.loc[:,'2020_01_COVID']+=df.loc[:,'2020_02_COVID']+df.loc[:,'2020_03_COVID']+df.loc[:,'2020_04_COVID']+df.loc[:,'2020_05_COVID']+df.loc[:,'2020_06_COVID']
        current_df.loc[:,'2020_01_non_COVID']+=df.loc[:,'2020_02_non_COVID']+df.loc[:,'2020_03_non_COVID']+df.loc[:,'2020_04_non_COVID']+df.loc[:,'2020_05_non_COVID']+df.loc[:,'2020_06_non_COVID']
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2, p, dof))


        current_df = df.loc[:, ["2016_01", "2017_01"]].copy()
        current_df.loc[:, '2016_01'] += df.loc[:, '2016_02'] + df.loc[:, '2016_03'] + df.loc[:, '2016_04'] + \
                                        df.loc[:, '2016_05'] + df.loc[:, '2016_06']
        current_df.loc[:, '2017_01'] += df.loc[:, '2017_02'] + df.loc[:, '2017_03'] + df.loc[:, '2017_04'] + \
                                        df.loc[:, '2017_05'] + df.loc[:, '2017_06']
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("2016/2017 chi2 {}, p {}".format(chi2, p, dof))

        current_df = df.loc[:, ["2017_01", "2018_01"]].copy()
        current_df.loc[:, '2017_01'] += df.loc[:, '2017_02'] + df.loc[:, '2017_03'] + df.loc[:, '2017_04'] + \
                                        df.loc[:, '2017_05'] + df.loc[:, '2017_06']
        current_df.loc[:, '2018_01'] += df.loc[:, '2018_02'] + df.loc[:, '2018_03'] + df.loc[:, '2018_04'] + \
                                        df.loc[:, '2018_05'] + df.loc[:, '2018_06']
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("2017/2018 chi2 {}, p {}".format(chi2, p, dof))

        current_df = df.loc[:, ["2018_01", "2019_01"]].copy()
        current_df.loc[:, '2018_01'] += df.loc[:, '2018_02'] + df.loc[:, '2018_03'] + df.loc[:, '2018_04'] + \
                                    df.loc[:, '2018_05'] + df.loc[:, '2018_06']
        current_df.loc[:, '2019_01'] += df.loc[:, '2019_02'] + df.loc[:, '2019_03'] + df.loc[:, '2019_04'] + \
                                    df.loc[:, '2019_05'] + df.loc[:, '2019_06']
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("2018/2019 chi2 {}, p {}".format(chi2, p, dof))

        current_df = df.loc[:, ["2019_01", "2020_01"]].copy()
        current_df.loc[:, '2019_01'] += df.loc[:, '2019_02'] + df.loc[:, '2019_03'] + df.loc[:, '2019_04'] + \
                                        df.loc[:, '2019_05'] + df.loc[:, '2019_06']
        current_df.loc[:, '2020_01'] += df.loc[:, '2020_02'] + df.loc[:, '2020_03'] + df.loc[:, '2020_04'] + \
                                        df.loc[:, '2020_05'] + df.loc[:, '2020_06']
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("2019/2020 chi2 {}, p {}".format(chi2, p, dof))

        current_df = df.loc[:, ["2019_01", "2020_01_non_COVID"]].copy()
        current_df.loc[:, '2019_01'] += df.loc[:, '2019_02'] + df.loc[:, '2019_03'] + df.loc[:, '2019_04'] + \
                                        df.loc[:, '2019_05'] + df.loc[:, '2019_06']
        current_df.loc[:,'2020_01_non_COVID']+=df.loc[:,'2020_02_non_COVID']+df.loc[:,'2020_03_non_COVID']+\
                                        df.loc[:,'2020_04_non_COVID']+df.loc[:,'2020_05_non_COVID']+df.loc[:,'2020_06_non_COVID']
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("2019/2020 Non COVID chi2 {}, p {}".format(chi2, p, dof))

        current_df = df.loc[:, ["2019_01", "2020_01_COVID"]].copy()
        current_df.loc[:, '2019_01'] += df.loc[:, '2019_02'] + df.loc[:, '2019_03'] + df.loc[:, '2019_04'] + \
                                        df.loc[:, '2019_05'] + df.loc[:, '2019_06']
        current_df.loc[:, '2020_01_COVID'] += df.loc[:, '2020_02_COVID'] + df.loc[:, '2020_03_COVID'] + \
                                                  df.loc[:, '2020_04_COVID'] + df.loc[:,'2020_05_COVID'] + df.loc[:,'2020_06_COVID']
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("2019/2020 COVID chi2 {}, p {}".format(chi2, p, dof))

        current_df = df.loc[:, ["2019_01", "2020_01_non_COVID"]].copy()
        current_df.loc[:, '2019_01'] += df.loc[:, '2019_02'] + df.loc[:, '2019_03'] + df.loc[:, '2019_04'] + \
                                        df.loc[:, '2019_05'] + df.loc[:, '2019_06']
        current_df.loc[:, '2020_01_non_COVID'] += df.loc[:, '2020_02_non_COVID'] + df.loc[:, '2020_03_non_COVID'] + \
                                              df.loc[:, '2020_04_non_COVID'] + df.loc[:, '2020_05_non_COVID'] + df.loc[:,
                                                                                                        '2020_06_non_COVID']
        current_df = current_df.T
        chi2, p, dof, expected = stats.chi2_contingency(current_df)
        print(current_df)
        print("2019/2020 non COVID chi2 {}, p {}".format(chi2, p, dof))

    def run_chi_single_vs_multi_authors(self, df):
        a_single=len(df.loc[df.loc[:, '2020_02_COVID'] == 1])
        a_multi=len(df.loc[df.loc[:, '2020_02_COVID'] > 1])
        b_single=len(df.loc[df.loc[:, '2020_02_non_COVID'] == 1])
        b_multi=len(df.loc[df.loc[:, '2020_02_non_COVID'] > 1])
        observed_array=np.array([[a_single, a_multi], [b_single, b_multi]])
        print(observed_array)
        chi2, p, dof, expected = stats.chi2_contingency(observed_array)
        print("Feb -COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2, p))

        a_single = len(df.loc[df.loc[:, '2020_03_COVID'] == 1])
        a_multi = len(df.loc[df.loc[:, '2020_03_COVID'] > 1])
        b_single = len(df.loc[df.loc[:, '2020_03_non_COVID'] == 1])
        b_multi = len(df.loc[df.loc[:, '2020_03_non_COVID'] > 1])
        observed_array = np.array([[a_single, a_multi], [b_single, b_multi]])
        print(observed_array)
        chi2, p, dof, expected = stats.chi2_contingency(observed_array)
        print("Mar -COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2, p))

        a_single = len(df.loc[df.loc[:, '2020_04_COVID'] == 1])
        a_multi = len(df.loc[df.loc[:, '2020_04_COVID'] > 1])
        b_single = len(df.loc[df.loc[:, '2020_04_non_COVID'] == 1])
        b_multi = len(df.loc[df.loc[:, '2020_04_non_COVID'] > 1])
        observed_array = np.array([[a_single, a_multi], [b_single, b_multi]])
        print(observed_array)
        chi2, p, dof, expected = stats.chi2_contingency(observed_array)
        print("Apr -COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2, p))

        a_single = len(df.loc[df.loc[:, '2020_05_COVID'] == 1])
        a_multi = len(df.loc[df.loc[:, '2020_05_COVID'] > 1])
        b_single = len(df.loc[df.loc[:, '2020_05_non_COVID'] == 1])
        b_multi = len(df.loc[df.loc[:, '2020_05_non_COVID'] > 1])
        observed_array = np.array([[a_single, a_multi], [b_single, b_multi]])
        print(observed_array)
        chi2, p, dof, expected = stats.chi2_contingency(observed_array)
        print("May -COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2, p))

        a_single = len(df.loc[df.loc[:, '2020_06_COVID'] == 1])
        a_multi = len(df.loc[df.loc[:, '2020_06_COVID'] > 1])
        b_single = len(df.loc[df.loc[:, '2020_06_non_COVID'] == 1])
        b_multi = len(df.loc[df.loc[:, '2020_06_non_COVID'] > 1])
        observed_array = np.array([[a_single, a_multi], [b_single, b_multi]])
        print(observed_array)
        chi2, p, dof, expected = stats.chi2_contingency(observed_array)
        print("June -COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2,p))

        a_single = len(df.loc[df.loc[:, '2020_01_COVID'] == 1])+len(df.loc[df.loc[:, '2020_02_COVID'] == 1])+len(df.loc[df.loc[:, '2020_03_COVID'] == 1])+\
                   len(df.loc[df.loc[:, '2020_04_COVID'] == 1])+len(df.loc[df.loc[:, '2020_05_COVID'] == 1])+len(df.loc[df.loc[:, '2020_06_COVID'] == 1])
        a_multi = len(df.loc[df.loc[:, '2020_01_COVID'] > 1])+len(df.loc[df.loc[:, '2020_02_COVID'] > 1])+len(df.loc[df.loc[:, '2020_03_COVID'] > 1])+ \
                  len(df.loc[df.loc[:, '2020_04_COVID'] > 1])+len(df.loc[df.loc[:, '2020_05_COVID'] > 1])+len(df.loc[df.loc[:, '2020_06_COVID'] > 1])
        b_single = len(df.loc[df.loc[:, '2020_01_non_COVID'] == 1])+len(df.loc[df.loc[:, '2020_02_non_COVID'] == 1])+len(df.loc[df.loc[:, '2020_03_non_COVID'] == 1])+ \
                   len(df.loc[df.loc[:, '2020_04_non_COVID'] == 1])+len(df.loc[df.loc[:, '2020_05_non_COVID'] == 1])+len(df.loc[df.loc[:, '2020_06_non_COVID'] == 1])
        b_multi = len(df.loc[df.loc[:, '2020_01_non_COVID'] > 1])+len(df.loc[df.loc[:, '2020_02_non_COVID'] > 1])+len(df.loc[df.loc[:, '2020_03_non_COVID'] > 1])+ \
                  len(df.loc[df.loc[:, '2020_04_non_COVID'] > 1])+len(df.loc[df.loc[:, '2020_05_non_COVID'] > 1])+len(df.loc[df.loc[:, '2020_06_non_COVID'] > 1])
        observed_array = np.array([[a_single, a_multi], [b_single, b_multi]])
        print(observed_array)
        chi2, p, dof, expected = stats.chi2_contingency(observed_array)
        print("COVID-19/Non COVID-19 chi2 {}, p {}".format(chi2, p, dof))


        a_single = len(df.loc[df.loc[:, '2016_01'] == 1])+len(df.loc[df.loc[:, '2016_02'] == 1])+len(df.loc[df.loc[:, '2016_03'] == 1])+ \
                   len(df.loc[df.loc[:, '2016_04'] == 1])+len(df.loc[df.loc[:, '2016_05'] == 1])+len(df.loc[df.loc[:, '2016_06'] == 1])
        a_multi = len(df.loc[df.loc[:, '2016_01'] > 1])+len(df.loc[df.loc[:, '2016_02'] > 1])+len(df.loc[df.loc[:, '2016_03'] > 1])+ \
                  len(df.loc[df.loc[:, '2016_04'] > 1])+len(df.loc[df.loc[:, '2016_05'] > 1])+len(df.loc[df.loc[:, '2016_06'] > 1])
        b_single = len(df.loc[df.loc[:, '2017_01'] == 1])+len(df.loc[df.loc[:, '2017_02'] == 1])+len(df.loc[df.loc[:, '2017_03'] == 1])+ \
                   len(df.loc[df.loc[:, '2017_04'] == 1])+len(df.loc[df.loc[:, '2017_05'] == 1])+len(df.loc[df.loc[:, '2017_06'] == 1])
        b_multi = len(df.loc[df.loc[:, '2017_01'] > 1])+len(df.loc[df.loc[:, '2017_02'] > 1])+len(df.loc[df.loc[:, '2017_03'] > 1])+ \
                  len(df.loc[df.loc[:, '2017_04'] > 1])+len(df.loc[df.loc[:, '2017_05'] > 1])+len(df.loc[df.loc[:, '2017_06'] > 1])
        observed_array = np.array([[a_single, a_multi], [b_single, b_multi]])
        print(observed_array)
        chi2, p, dof, expected = stats.chi2_contingency(observed_array)
        print("2016/2017 chi2 {}, p {}".format(chi2, p, dof))


        a_single = len(df.loc[df.loc[:, '2017_01'] == 1])+len(df.loc[df.loc[:, '2017_02'] == 1])+len(df.loc[df.loc[:, '2017_03'] == 1])+ \
                   len(df.loc[df.loc[:, '2017_04'] == 1])+len(df.loc[df.loc[:, '2017_05'] == 1])+len(df.loc[df.loc[:, '2017_06'] == 1])
        a_multi = len(df.loc[df.loc[:, '2017_01'] > 1])+len(df.loc[df.loc[:, '2017_02'] > 1])+len(df.loc[df.loc[:, '2017_03'] > 1])+ \
                  len(df.loc[df.loc[:, '2017_04'] > 1])+len(df.loc[df.loc[:, '2017_05'] > 1])+len(df.loc[df.loc[:, '2017_06'] > 1])
        b_single = len(df.loc[df.loc[:, '2018_01'] == 1])+len(df.loc[df.loc[:, '2018_02'] == 1])+len(df.loc[df.loc[:, '2018_03'] == 1])+ \
                   len(df.loc[df.loc[:, '2018_04'] == 1])+len(df.loc[df.loc[:, '2018_05'] == 1])+len(df.loc[df.loc[:, '2018_06'] == 1])
        b_multi = len(df.loc[df.loc[:, '2018_01'] > 1])+len(df.loc[df.loc[:, '2018_02'] > 1])+len(df.loc[df.loc[:, '2018_03'] > 1])+ \
                  len(df.loc[df.loc[:, '2018_04'] > 1])+len(df.loc[df.loc[:, '2018_05'] > 1])+len(df.loc[df.loc[:, '2018_06'] > 1])
        observed_array = np.array([[a_single, a_multi], [b_single, b_multi]])
        print(observed_array)
        chi2, p, dof, expected = stats.chi2_contingency(observed_array)
        print("2017/2018 chi2 {}, p {}".format(chi2, p, dof))

        a_single = len(df.loc[df.loc[:, '2018_01'] == 1]) + len(df.loc[df.loc[:, '2018_02'] == 1]) + len(df.loc[df.loc[:, '2018_03'] == 1]) + \
                   len(df.loc[df.loc[:, '2018_04'] == 1]) + len(df.loc[df.loc[:, '2018_05'] == 1]) + len(df.loc[df.loc[:, '2018_06'] == 1])
        a_multi = len(df.loc[df.loc[:, '2018_01'] > 1]) + len(df.loc[df.loc[:, '2018_02'] > 1]) + len(df.loc[df.loc[:, '2018_03'] > 1]) + \
                  len(df.loc[df.loc[:, '2018_04'] > 1]) + len(df.loc[df.loc[:, '2018_05'] > 1]) + len(df.loc[df.loc[:, '2018_06'] > 1])
        b_single = len(df.loc[df.loc[:, '2019_01'] == 1]) + len(df.loc[df.loc[:, '2019_02'] == 1]) + len(df.loc[df.loc[:, '2019_03'] == 1]) + \
                   len(df.loc[df.loc[:, '2019_04'] == 1]) + len(df.loc[df.loc[:, '2019_05'] == 1]) + len(df.loc[df.loc[:, '2019_06'] == 1])
        b_multi = len(df.loc[df.loc[:, '2019_01'] > 1]) + len(df.loc[df.loc[:, '2019_02'] > 1]) + len(df.loc[df.loc[:, '2019_03'] > 1]) + \
                  len(df.loc[df.loc[:, '2019_04'] > 1]) + len(df.loc[df.loc[:, '2019_05'] > 1]) + len(df.loc[df.loc[:, '2019_06'] > 1])
        observed_array = np.array([[a_single, a_multi], [b_single, b_multi]])
        print(observed_array)
        chi2, p, dof, expected = stats.chi2_contingency(observed_array)
        print("2018/2019 chi2 {}, p {}".format(chi2, p, dof))

        a_single = len(df.loc[df.loc[:, '2019_01'] == 1]) + len(df.loc[df.loc[:, '2019_02'] == 1]) + len(df.loc[df.loc[:, '2019_03'] == 1]) + \
                   len(df.loc[df.loc[:, '2019_04'] == 1]) + len(df.loc[df.loc[:, '2019_05'] == 1]) + len(df.loc[df.loc[:, '2019_06'] == 1])
        a_multi = len(df.loc[df.loc[:, '2019_01'] > 1]) + len(df.loc[df.loc[:, '2019_02'] > 1]) + len(df.loc[df.loc[:, '2019_03'] > 1]) + \
                  len(df.loc[df.loc[:, '2019_04'] > 1]) + len(df.loc[df.loc[:, '2019_05'] > 1]) + len(df.loc[df.loc[:, '2019_06'] > 1])
        b_single = len(df.loc[df.loc[:, '2020_01'] == 1]) + len(df.loc[df.loc[:, '2020_02'] == 1]) + len(df.loc[df.loc[:, '2020_03'] == 1]) + \
                   len(df.loc[df.loc[:, '2020_04'] == 1]) + len(df.loc[df.loc[:, '2020_05'] == 1]) + len(df.loc[df.loc[:, '2020_06'] == 1])
        b_multi = len(df.loc[df.loc[:, '2020_01'] > 1]) + len(df.loc[df.loc[:, '2020_02'] > 1]) + len(df.loc[df.loc[:, '2020_03'] > 1]) + \
                  len(df.loc[df.loc[:, '2020_04'] > 1]) + len(df.loc[df.loc[:, '2020_05'] > 1]) + len(df.loc[df.loc[:, '2020_06'] > 1])
        observed_array = np.array([[a_single, a_multi], [b_single, b_multi]])
        print(observed_array)
        chi2, p, dof, expected = stats.chi2_contingency(observed_array)
        print("2019/2020 chi2 {}, p {}".format(chi2, p, dof))

        a_single = len(df.loc[df.loc[:, '2019_01'] == 1]) + len(df.loc[df.loc[:, '2019_02'] == 1]) + len(df.loc[df.loc[:, '2019_03'] == 1]) + \
                   len(df.loc[df.loc[:, '2019_04'] == 1]) + len(df.loc[df.loc[:, '2019_05'] == 1]) + len(df.loc[df.loc[:, '2019_06'] == 1])
        a_multi = len(df.loc[df.loc[:, '2019_01'] > 1]) + len(df.loc[df.loc[:, '2019_02'] > 1]) + len(df.loc[df.loc[:, '2019_03'] > 1]) + \
                  len(df.loc[df.loc[:, '2019_04'] > 1]) + len(df.loc[df.loc[:, '2019_05'] > 1]) + len(df.loc[df.loc[:, '2019_06'] > 1])
        b_single = len(df.loc[df.loc[:, '2020_01_non_COVID'] == 1]) + len(df.loc[df.loc[:, '2020_02_non_COVID'] == 1]) + len(df.loc[df.loc[:, '2020_03_non_COVID'] == 1]) + \
                   len(df.loc[df.loc[:, '2020_04_non_COVID'] == 1]) + len(df.loc[df.loc[:, '2020_05_non_COVID'] == 1]) + len(df.loc[df.loc[:, '2020_06_non_COVID'] == 1])
        b_multi = len(df.loc[df.loc[:, '2020_01_non_COVID'] > 1]) + len(df.loc[df.loc[:, '2020_02_non_COVID'] > 1]) + len(df.loc[df.loc[:, '2020_03_non_COVID'] > 1]) + \
                  len(df.loc[df.loc[:, '2020_04_non_COVID'] > 1]) + len(df.loc[df.loc[:, '2020_05_non_COVID'] > 1]) + len(df.loc[df.loc[:, '2020_06_non_COVID'] > 1])
        observed_array = np.array([[a_single, a_multi], [b_single, b_multi]])
        print(observed_array)
        chi2, p, dof, expected = stats.chi2_contingency(observed_array)
        print("2019/2020 Non Covid chi2 {}, p {}".format(chi2, p, dof))


    def get_pearsons_correlation_covid_counts(self,df):
        # df=df.drop(['New Scientist', 'Medical Hypotheses','Diabetes & Metabolic Syndrome: Clinical Research & Reviews'])
        # df=df.drop(['New Scientist'])

        df['COVID-19 papers'] = df['2020-01_COVID'] + df['2020-02_COVID']+df['2020-03_COVID'] + df['2020-04_COVID']+df['2020-05_COVID'] + df['2020-06_COVID']
        df['sum']=df['2020-01'] + df['2020-02']+df['2020-03'] + df['2020-04']+df['2020-05'] + df['2020-06']
        df_calc=df[['COVID-19 papers','SJR', 'sum']].copy().T
        df_calc.loc['sum', 'Brazilian Journal of Anesthesiology'] = 62
        df_calc.loc['COVID-19 papers', 'Brazilian Journal of Anesthesiology'] = 7
        df_calc.loc['SJR', 'Brazilian Journal of Anesthesiology'] = 0.11
        df_calc.loc['sum', 'Journal of Dental Sciences'] = 157
        df_calc.loc['COVID-19 papers', 'Journal of Dental Sciences'] = 11
        df_calc.loc['SJR', 'Journal of Dental Sciences'] = 0.23
        df_calc.loc['sum', 'Pulmonology'] = 112
        df_calc.loc['COVID-19 papers', 'Pulmonology'] = 18
        df_calc.loc['SJR', 'Pulmonology'] = 0.56
        df_calc.loc['sum', 'Visual Journal of Emergency Medicine'] = 148
        df_calc.loc['COVID-19 papers', 'Visual Journal of Emergency Medicine'] = 7
        df_calc.loc['SJR', 'Visual Journal of Emergency Medicine'] = 0.11
        df_calc.loc['sum', 'Indian Journal of Tuberculosis'] = 87
        df_calc.loc['COVID-19 papers', 'Indian Journal of Tuberculosis'] = 3
        df_calc.loc['SJR', 'Indian Journal of Tuberculosis'] = 0.294
        df_calc.loc['sum', 'Medical Mycology Case Reports'] = 52
        df_calc.loc['COVID-19 papers', 'Medical Mycology Case Reports'] = 3
        df_calc.loc['SJR', 'Medical Mycology Case Reports'] = 0.404
        df_calc.loc['sum', 'Comparative Immunology, Microbiology and Infectious Diseases'] = 91
        df_calc.loc['COVID-19 papers', 'Comparative Immunology, Microbiology and Infectious Diseases'] = 1
        df_calc.loc['SJR', 'Comparative Immunology, Microbiology and Infectious Diseases'] = 0.607
        #
        df_calc.loc['sum', 'The Lancet'] = 943
        df_calc.loc['COVID-19 papers', 'The Lancet'] = 344
        df_calc.loc['SJR', 'The Lancet'] = 14.55
        df_calc.loc['sum', 'Journal of the American Academy of Dermatology'] = 1252
        df_calc.loc['COVID-19 papers', 'Journal of the American Academy of Dermatology'] = 148
        df_calc.loc['SJR', 'Journal of the American Academy of Dermatology'] = 1.81
        df_calc.loc['sum', 'Science of The Total Environment'] = 6044
        df_calc.loc['COVID-19 papers', 'Science of The Total Environment'] = 156
        df_calc.loc['SJR', 'Science of The Total Environment'] = 1.66

        # df_calc.loc['COVID-19 papers', 'cell'] = 38
        # df_calc.loc['SJR', 'cell'] = 24.698
        # df_calc.loc['COVID-19 papers', 'The Lancet oncology'] = 40
        # df_calc.loc['SJR', 'The Lancet oncology'] = 15.65
        # df_calc.loc['COVID-19 papers', 'Immunity'] = 20
        # df_calc.loc['SJR', 'Immunity'] = 11.977
        # df_calc.loc['COVID-19 papers', 'journal of clinical virology'] = 138
        # df_calc.loc['SJR', 'journal of clinical virology'] = 1.5
        # df_calc.loc['COVID-19 papers', 'world neurosurgery'] = 104
        # df_calc.loc['SJR', 'world neurosurgery'] = 0.73
        # df_calc.loc['COVID-19 papers', 'american journal of infection control'] = 68
        # df_calc.loc['SJR', 'american journal of infection control'] = 0.99
        # df_calc.loc['COVID-19 papers', 'gastroenterology'] = 109
        # df_calc.loc['SJR', 'gastroenterology'] = 6.85
        # df_calc.loc['COVID-19 papers', 'diabetes research and clinical practice'] = 54
        # df_calc.loc['SJR', 'diabetes research and clinical practice'] = 1.4

        r,pValue=stats.pearsonr(df_calc.T['COVID-19 papers'], df_calc.T['SJR'])
        # r,pValue=stats.pearsonr(df['COVID-19 papers'], df['SJR'])
        print('COVID total Pearsons r {}. pValue {}'.format(r, pValue))

        r, pValue = stats.pearsonr(df['COVID-19 papers'], df['SJR'])
        # r,pValue=stats.pearsonr(df['COVID-19 papers'], df['SJR'])
        print('COVID selected journals Pearsons r {}. pValue {}'.format(r, pValue))

        r, pValue = stats.pearsonr(df['2020-06_COVID'], df['SJR'])
        print('COVID June Pearsons r {}. pValue {}'.format(r, pValue))
        df['sum']=df['2020-01'] + df['2020-02']+df['2020-03'] + df['2020-04']+df['2020-05'] + df['2020-06']
        df['percent']=100*df['COVID-19 papers']/df['sum']
        r, pValue = stats.pearsonr(df['percent'], df['SJR'])
        print('COVID out of total per journal percent Pearsons r {}. pValue {}'.format(r, pValue))

        sum_subset=df['COVID-19 papers'].sum()
        df['percent_COVID'] = 100 * df['COVID-19 papers'] / sum_subset
        r, pValue = stats.pearsonr(df['percent_COVID'], df['SJR'])
        print('COVID subset percent Pearsons r {}. pValue {}'.format(r, pValue))

        # df_calc.T.plot(x='SJR', y='COVID-19 papers', kind='scatter', color='G', logx=True, xlim=(0,20))
        # df_calc.T.plot.scatter(x='SJR', y='COVID-19 papers')
        # df.plot(x='SJR', y='COVID-19 papers', kind='scatter', color='R', logx=True)

        df_calc=df_calc.T
        df_calc['percent'] = 100 * df_calc['COVID-19 papers'] / df_calc['sum']
        r, pValue = stats.pearsonr(df_calc['percent'], df_calc['SJR'])
        print('COVID out of total per journal full list percent Pearsons r {}. pValue {}'.format(r, pValue))

        sum=df_calc['COVID-19 papers'].sum()
        df_calc['percent_COVID']=100*df_calc['COVID-19 papers']/sum

        df_calc=df_calc.sort_values(by='SJR', ascending=False)
        print(df_calc)
        r, pValue = stats.pearsonr(df_calc['percent_COVID'], df_calc['SJR'])
        print('COVID percent full list Pearsons r {}. pValue {}'.format(r, pValue))

        # utils.write_to_csv(df_calc,'counts_by_SJR_COVID.csv', index=True)
        tex = df_calc.to_latex(index=True)
        text_file = open("scienceDirectData/counts_by_SJR_COVID.tex", "w")
        text_file.write(tex)
        text_file.close()



    def get_pearsons_correlation_covid_acc_time(self,df):
        # df = df.drop(
        #     ['Medical Hypotheses', 'Diabetes & Metabolic Syndrome: Clinical Research & Reviews'])
        cov_means = df.T.loc[['2020-01_COVID',
                                            '2020-02_COVID',
                                            '2020-03_COVID',
                                            '2020-04_COVID',
                                            '2020-05_COVID',
                                            '2020-06_COVID',]]
        cov_means_mean= cov_means.mean()
        # df['mean_COVID'] = (df['2020-01_COVID'] + df['2020-02_COVID'] + df['2020-03_COVID'] + df['2020-04_COVID'] + df[
        #     '2020-05_COVID'] + df['2020-06_COVID'])/6
        r, pValue = stats.pearsonr(cov_means_mean.values, df['SJR'].values)
        print('COVID total Pearsons r {}. pValue {}'.format(r, pValue))

        r, pValue = stats.pearsonr(df['2020-06_COVID'], df['SJR'])
        print('COVID June Pearsons r {}. pValue {}'.format(r, pValue))





    def get_and_store_acc_df(self):
        acc_df=sdd.get_all_acceptance_time_covid_non_covid(utils)
        utils.save_obj(acc_df, 'acceptance_df')

    def get_and_store_countries_collab_df(self):
        countries_dict=utils.load_obj("countries_dict")
        collab_countries_df=sdd.get_all_countries_collab(countries_dict)
        utils.save_obj(collab_countries_df, 'collab_countries_df')

    def get_and_store_countries_collab_by_papers_df(self):
        countries_dict=utils.load_obj("countries_dict")
        collab_countries_by_num_papers_df=sdd.get_all_countries_collab_num_papers(countries_dict)
        utils.save_obj(collab_countries_by_num_papers_df, 'collab_countries_by_num_papers_df')

    def get_and_store_single_vs_multi_country(self):
        countries_dict = utils.load_obj("countries_dict")
        single_vs_multi_country_df = sdd.get_all_countries_single_vs_multi(countries_dict)
        utils.save_obj(single_vs_multi_country_df, 'single_vs_multi_country_df')

    def get_and_store_authors(self):
        authors_df=sdd.get_all_authors(utils)
        utils.save_obj(authors_df, 'authors_df')




if __name__ == '__main__':
    utils = Utils(path=path)
    sdd=SciencedirectData()

    st = StatsTests()


    st.run_acc_ttests()
    # st.get_and_store_countries_collab_df()
    # df = utils.load_obj('collab_countries_df')
    # df.drop(['2016', '2017', '2018', '2019', '2020', '2020_COVID', '2020_non_COVID'], axis=1)
    # df = df.apply(pd.to_numeric, errors='coerce')
    # df = df.replace(0, np.NaN)
    # df = df.T
    # st.run_collab_countries_ttests(df)


    # df = utils.load_obj('collab_countries_by_num_papers_df')
    # df = df.replace(0, np.NaN)
    # df = df.T
    # st.run_collab_countries_by_num_papers_ttests(df)

    # st.get_and_store_single_vs_multi_country()
    df = utils.load_obj('single_vs_multi_country_df')
    # df = df.replace(0, np.NaN)
    # st.run_chi_single_vs_multi_country(df)

    # st.get_and_store_authors()
    # df = utils.load_obj('authors_df')
    # df = df.replace(0, np.NaN)get_pearsons_correlation_covid_counts
    # df = df.T
    # st.run_chi_single_vs_multi_authors(df)

    df=utils.load_csv_data_to_df("journals_counts_mean.csv")
    df,df_COVID=sdd.create_counts_table(df, utils)
    st.get_pearsons_correlation_covid_counts(df_COVID)
    #
    # df=utils.load_csv_data_to_df("journals_counts_mean.csv")
    # df=sdd.get_acc_time_journals(df, utils)
    # st.get_pearsons_correlation_covid_acc_time(df)



