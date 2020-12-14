
import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap
from utils import Utils
import math
from tqdm import tqdm_notebook
# %matplotlib inline
from parse_scopus_search import ParseScopusSearch
import csv
import os

class Visualization():
    def show_journals(self, df):
        df.journal.value_counts()[0:10].plot(kind='bar')
        plt.grid()
        plt.show()

    def show_sources(self,df):
        df.source_x.value_counts().plot(kind='bar')
        plt.show()

    def show_year(self,df):
        df.publish_time.value_counts()[0:10].plot(kind='bar')
        plt.show()

    def plt_by_month(self, medrxiv_data, biorxiv_data=None, arxiv_data=None, scopus_data=None,world_infectious_statistics=None ):
        t=["Oct","Nov","Dec","jan","Feb","Mar","Apr"]

        # plt.subplot(4,1,1)
        plt.plot(t, medrxiv_data[0], color='blue', label='medrxiv')
        plt.plot(t, biorxiv_data[0], color='magenta', label='biorxiv')
        plt.plot(t, arxiv_data[0], color='green', label='arxiv;q-bio')
        plt.plot(t, scopus_data[0], color='black', label='scopus')
        plt.plot(t, world_infectious_statistics[0], color='red', label='confirmed virus infections worldwide, factor of 1/1000 ')
        plt.ylabel('# of published papers')
        plt.xlabel('Months in 2019, 2020')
        plt.title('number of COVID-19 related publications')
        plt.legend()
        plt.grid(True)
        plt.show()

        plt.close()

        plt.plot(t, medrxiv_data[1], color='blue', linestyle="--", label='medrxiv')
        plt.plot(t, biorxiv_data[1], color='magenta', linestyle="--", label='biorxiv')
        plt.plot(t, arxiv_data[1], color='green', linestyle="--", label='arxiv;q-bio')
        plt.ylabel('# of published papers')
        plt.xlabel('Months in 2019, 2020')
        plt.title('number of total publications')
        plt.legend()
        plt.show()

        plt.close()

        # plt.subplot(2,1,2)
        plt.plot(t,  [x * 100 for x in medrxiv_data[2]], color='blue', linestyle="dotted", label='medrxiv')
        plt.plot(t, [x * 100 for x in biorxiv_data[2]], color='magenta', linestyle="dotted", label='biorxiv')
        plt.plot(t, [x * 100 for x in arxiv_data[2]], color='green', linestyle="dotted", label='arxiv;q-bio')

        plt.xlabel('Months in 2019, 2020')
        plt.ylabel('covid-19 rate')
        plt.title('rate of COVID-19 related publications, shown in %')
        plt.legend()
        plt.grid(True)
        plt.show()


    def plt_journals_counts(self,df, title):
        f = plt.figure()
        ax = f.gca()
        plt.title(title, color='black')
        df.plot(ax=ax)
        # ax.set_xticks(df.index)
        ax.set_xticklabels(df.index, rotation=45)
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.show()

    def plt_journals_sums_old(self,COV, Totals, title):


        fig, axs = plt.subplots(1,5,sharey=True, facecolor='w')
        ax1=axs[0]
        ax2=axs[1]
        ax3=axs[2]
        ax4=axs[3]
        ax5=axs[4]
        ax6 = ax5.twinx()
        ax1.set_ylim(0,2000)
        ax6.set_ylim(0,800)
        ax1.set_xlim(-1,3)
        ax2.set_xlim(4,7)
        ax3.set_xlim(8,11)
        ax4.set_xlim(12,15)
        ax5.set_xlim(15.75,20)
        ax1.plot(*zip(*Totals.items()), color='blue', label='Total publications')
        ax2.plot(*zip(*Totals.items()), color='blue', label='Total publications')
        ax3.plot(*zip(*Totals.items()), color='blue', label='Total publications')
        ax4.plot(*zip(*Totals.items()), color='blue', label='Total publications')
        ax5.plot(*zip(*Totals.items()), color='blue', label='Total publications')


        ax1.spines['right'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.spines['left'].set_visible(False)
        ax5.spines['right'].set_visible(False)
        ax5.spines['left'].set_visible(False)
        ax6.spines['left'].set_visible(False)
        ax1.yaxis.tick_left()

        ax1.tick_params(axis='y',color='white')
        ax2.tick_params(axis='y',color='white')
        ax3.tick_params(axis='y',color='white')
        ax4.tick_params(axis='y',color='white')
        ax5.tick_params(axis='y',color='white')

        ax6.yaxis.tick_right()


        ax6.bar(*zip(*COV.items()), width=0.1, color='green', label='COVID-19 publications')

        d = .005  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (-d, +d), **kwargs)
        ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)


        kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
        ax3.plot((-d, +d), (-d, +d), **kwargs)
        ax3.plot((1 - d, 1 + d), (-d, +d), **kwargs)


        kwargs.update(transform=ax4.transAxes)  # switch to the bottom axes
        ax4.plot((-d, +d), (-d, +d), **kwargs)
        ax4.plot((1 - d, 1 + d), (-d, +d), **kwargs)


        kwargs.update(transform=ax5.transAxes)  # switch to the bottom axes
        ax5.plot((-d, +d), (-d, +d), **kwargs)

        # ax=plt.bar(*zip(*COV.items()), width=0.1, color='green', label='COVID-19 publications')
        # plt.plot(secondary_y=True, *zip(*Totals.items()), label='Total publications', ax=ax)

        # ax1.legend()
        # ax6.legend()
        # plt.xticks(range(len(COV)), list(COV.keys()), rotation=45)


        ax3.set_xlabel('First 4 months in 2016-2020')
        ax1.set_ylabel('Total publication counts', color='b')
        ax6.set_ylabel('COVID-19 related publication counts', color='g')

        # ax.set_xticks(df.index)
        ax1.set_xticklabels(list(COV.keys()), rotation=70)
        ax2.set_xticklabels(list(COV.keys()), rotation=70)
        ax3.set_xticklabels(list(COV.keys()), rotation=70)
        ax4.set_xticklabels(list(COV.keys()), rotation=70)
        ax5.set_xticklabels(list(COV.keys()), rotation=70)
        ax6.set_xticklabels(list(COV.keys()), rotation=70)
        ax3.set_title(title, color='black')

        plt.show()


    def plt_journals_sums(self,df, cov_df, non_cov_df,title):


        fig, axs = plt.subplots(1,5,sharey=True, facecolor='w')
        ax1=axs[0]
        ax2=axs[1]
        ax3=axs[2]
        ax4=axs[3]
        ax5=axs[4]
        ax6 = ax5.twinx()
        ax7 = ax5.twiny()
        ax1.set_ylim(0,900)
        ax6.set_ylim(0,600)
        ax1.set_xlim(-1,5)
        ax2.set_xlim(6,11)
        ax3.set_xlim(12,17)
        ax4.set_xlim(18,23)
        ax5.set_xlim(24,29.4)

        br1 = range(len(df['Total']))
        ax1.plot(br1, list(df['Total'].values), color='magenta',linewidth=2, label='Total publications')
        ax2.plot(br1, list(df['Total'].values), color='magenta', linewidth=2)
        ax3.plot(br1, list(df['Total'].values), color='magenta', linewidth=2)
        ax4.plot(br1, list(df['Total'].values), color='magenta', linewidth=2)
        ax5.plot(br1, list(df['Total'].values), color='magenta', linewidth=2)


        ax1.spines['right'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.spines['left'].set_visible(False)
        ax5.spines['right'].set_visible(False)
        ax5.spines['left'].set_visible(False)
        ax6.spines['left'].set_visible(False)
        ax7.spines['left'].set_visible(False)


        ax1.yaxis.tick_left()

        ax1.tick_params(axis='y',color='white', labelsize='16')
        ax2.tick_params(axis='y',color='white')
        ax3.tick_params(axis='y',color='white')
        ax4.tick_params(axis='y',color='white')
        ax5.tick_params(axis='y',color='white')

        ax6.yaxis.tick_right()
        ax6.tick_params(axis='y', color='white', labelsize='16')

        ax6.bar(br1, list(cov_df['Total'].values), width=0.1, color='green', label='COVID-19 publications')
        ax7.plot(br1, list(non_cov_df['Total'].values), color='blue', linewidth=2, label='Non COVID-19 publications')

        d = .005  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (-d, +d), **kwargs)
        ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)


        kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
        ax3.plot((-d, +d), (-d, +d), **kwargs)
        ax3.plot((1 - d, 1 + d), (-d, +d), **kwargs)


        kwargs.update(transform=ax4.transAxes)  # switch to the bottom axes
        ax4.plot((-d, +d), (-d, +d), **kwargs)
        ax4.plot((1 - d, 1 + d), (-d, +d), **kwargs)


        kwargs.update(transform=ax5.transAxes)  # switch to the bottom axes
        ax5.plot((-d, +d), (-d, +d), **kwargs)

        # ax=plt.bar(*zip(*COV.items()), width=0.1, color='green', label='COVID-19 publications')
        # plt.plot(secondary_y=True, *zip(*Totals.items()), label='Total publications', ax=ax)

        # ax1.legend()
        # ax6.legend()
        # plt.xticks(range(len(COV)), list(COV.keys()), rotation=45)


        ax3.set_xlabel('First six months in 2016-2020', fontdict={'fontsize': 16})
        ax1.set_ylabel('Total publication counts', fontdict={'fontsize': 16}, color='magenta')
        ax6.set_ylabel('COVID-19 and non COVID-19 related publication counts', fontdict={'fontsize': 16}, color='g')


        # ax.set_xticks(df.index)
        ax1.set_xticklabels(list(df.T.keys()), rotation=70, fontdict={'fontsize': 16})
        ax2.set_xticklabels(list(df.T.keys()), rotation=70, fontdict={'fontsize': 16})
        ax3.set_xticklabels(list(df.T.keys()), rotation=70, fontdict={'fontsize': 16})
        ax4.set_xticklabels(list(df.T.keys()), rotation=70, fontdict={'fontsize': 16})
        ax5.set_xticklabels(list(df.T.keys()), rotation=70, fontdict={'fontsize': 16})
        ax6.set_xticklabels(list(df.T.keys()), rotation=70, fontdict={'fontsize': 16})
        ax7.set_xticklabels([])
        ax3.set_title(title, color='black', fontdict={'fontsize': 18})
        fig.legend(loc=(0.08, 0.70), borderaxespad=0.1, fontsize=14)
        plt.show()

    def plt_rxivs_cov_and_totals_old(self,df, title):
        barWidth = 0.1
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax2.set_ylim(0,2000)

        br1=range(len(df['medrxiv_cov']))
        br2=[x + barWidth for x in br1]
        br3 = [x + barWidth for x in br2]
        # br4 = [x + barWidth for x in br3]
        ax1.plot(br1,list(df['medrxiv'].values), color='blue', linestyle='dashed',label='medriv total')
        # ax2.bar(*zip(*df['medrxiv_cov'].items()), width=0.05, color='blue', label='COVID-19 publications medrxiv')
        ax2.bar(br1, list(df['medrxiv_cov'].values), width=0.1, color='blue', label='COVID-19 publications medrxiv')

        ax1.plot(br1,list(df['biorxiv'].values), color='green', linestyle='dashed', label='biorxiv total')
        ax2.bar(br2, list(df['biorxiv_cov'].values), width=0.1, color='green', label='COVID-19 publications biorxiv')

        ax1.plot(br1,list(df['arxiv'].values), color='magenta', linestyle='dashed', label='arxiv (q-bio) total')
        ax2.bar(br3,list(df['arxiv_cov'].values), width=0.1, color='magenta', label='COVID-19 publications arxiv')

        ax1.plot(br1,list(df['world_stats'].values), color='red', label='confirmed virus infections worldwide, factor of 1/1000 ')

        ax1.set_xticklabels(list(df.T.keys()), rotation=70)
        ax1.set_xlabel('First 4 months in 2016-2020')
        ax1.set_ylabel('Total publication counts', color='b')
        ax2.set_ylabel('COVID-19 related publication counts', color='g')

        ax1.legend(loc=(0.01, 0.75))
        ax2.legend(loc=(0.01, 0.55))
        plt.title(title, color='black')

        plt.show()

    def plt_rxivs_cov_and_totals(self,df, title):
        barWidth = 0.1
        fig, axs = plt.subplots(1, 5, sharey=True, facecolor='w')
        ax1 = axs[0]
        ax2 = axs[1]
        ax3 = axs[2]
        ax4 = axs[3]
        ax5 = axs[4]
        ax6 = ax5.twinx()
        # fig, ax1 = plt.subplots()
        # ax2 = ax1.twinx()
        ax1.set_ylim(0,5000)
        ax6.set_ylim(0,2500)

        ax1.set_xlim(0, 5)
        ax2.set_xlim(6, 11)
        ax3.set_xlim(12, 17)
        ax4.set_xlim(18, 23)
        ax5.set_xlim(24, 29.6)

        br1=range(len(df['medrxiv_cov']))
        br2=[x + barWidth for x in br1]
        br3 = [x + barWidth for x in br2]
        br4 = [x + barWidth for x in br3]
        ax1.plot(br1,list(df['medrxiv'].values), color='blue', linestyle='dashed',linewidth=2, label='medrxiv total')
        ax2.plot(br1, list(df['medrxiv'].values), color='blue', linestyle='dashed', linewidth=2)
        ax3.plot(br1, list(df['medrxiv'].values), color='blue', linestyle='dashed', linewidth=2)
        ax4.plot(br1, list(df['medrxiv'].values), color='blue', linestyle='dashed', linewidth=2)
        ax5.plot(br1, list(df['medrxiv'].values), color='blue', linestyle='dashed', linewidth=2)

        ax6.bar(br3, list(df['medrxiv_cov'].values), width=0.1, color='blue', label='COVID-19 publications medrxiv')

        ax1.plot(br1,list(df['biorxiv'].values), color='green', linestyle='dashed', linewidth=2, label='biorxiv total')
        ax2.plot(br1, list(df['biorxiv'].values), color='green', linestyle='dashed', linewidth=2)
        ax3.plot(br1, list(df['biorxiv'].values), color='green', linestyle='dashed', linewidth=2)
        ax4.plot(br1, list(df['biorxiv'].values), color='green', linestyle='dashed', linewidth=2)
        ax5.plot(br1, list(df['biorxiv'].values), color='green', linestyle='dashed', linewidth=2)

        ax6.bar(br2, list(df['biorxiv_cov'].values), width=0.1, color='green', label='COVID-19 publications biorxiv')

        ax1.plot(br1,list(df['arxiv'].values), color='magenta', linestyle='dashed', linewidth=2, label='arxiv (q-bio) total')
        ax2.plot(br1,list(df['arxiv'].values), color='magenta', linestyle='dashed', linewidth=2)
        ax3.plot(br1,list(df['arxiv'].values), color='magenta', linestyle='dashed', linewidth=2)
        ax4.plot(br1,list(df['arxiv'].values), color='magenta', linestyle='dashed', linewidth=2)
        ax5.plot(br1,list(df['arxiv'].values), color='magenta', linestyle='dashed', linewidth=2)

        ax6.bar(br1,list(df['arxiv_cov'].values), width=0.1, color='magenta', label='COVID-19 publications arxiv')

        ax6.bar(br4,list(df['world_stats'].values), width=0.1, color='red', label='confirmed virus infections worldwide, factor of 1/5000 ')

        ax1.spines['right'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.spines['left'].set_visible(False)
        ax5.spines['right'].set_visible(False)
        ax5.spines['left'].set_visible(False)
        ax6.spines['left'].set_visible(False)
        ax1.yaxis.tick_left()

        ax1.tick_params(axis='y', color='white', labelsize='16')
        ax2.tick_params(axis='y', color='white')
        ax3.tick_params(axis='y', color='white')
        ax4.tick_params(axis='y', color='white')
        ax5.tick_params(axis='y', color='white')


        ax6.yaxis.tick_right()
        ax6.tick_params(axis='y', color='white', labelsize='16')

        ax1.set_xticklabels(list(df.T.keys()[0:6]), rotation=70, fontdict={'fontsize': 16})
        ax2.set_xticklabels(list(df.T.keys()[6:12]), rotation=70, fontdict={'fontsize': 16})
        ax3.set_xticklabels(list(df.T.keys()[12:18]), rotation=70,fontdict={'fontsize': 16})
        ax4.set_xticklabels(list(df.T.keys()[18:24]), rotation=70, fontdict={'fontsize': 16})
        ax5.set_xticklabels(list(df.T.keys()[24:30]), rotation=70, fontdict={'fontsize': 16})
        ax6.set_xticklabels(list(df.T.keys()[24:30]), rotation=70, fontdict={'fontsize': 16})

        d = .005  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (-d, +d), **kwargs)
        ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)

        kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
        ax3.plot((-d, +d), (-d, +d), **kwargs)
        ax3.plot((1 - d, 1 + d), (-d, +d), **kwargs)

        kwargs.update(transform=ax4.transAxes)  # switch to the bottom axes
        ax4.plot((-d, +d), (-d, +d), **kwargs)
        ax4.plot((1 - d, 1 + d), (-d, +d), **kwargs)

        kwargs.update(transform=ax5.transAxes)  # switch to the bottom axes
        ax5.plot((-d, +d), (-d, +d), **kwargs)

        ax3.set_xlabel('First six months in 2016-2020', fontdict={'fontsize': 16})

        ax1.set_ylabel('Total publication counts', fontdict={'fontsize': 16})
        ax6.set_ylabel('COVID-19 related publication counts', fontdict={'fontsize': 16})

        fig.legend(loc=(0.07, 0.65), borderaxespad=0.1, fontsize=14)
        ax3.set_title(title, color='black', fontdict={'fontsize': 18})

        plt.show()

    def plt_acc_time_avgs_history(self,df, df_sem, title):
        barWidth = 0.1
        fig, axs = plt.subplots(1, 5, sharey=True, facecolor='w')
        ax1 = axs[0]
        ax2 = axs[1]
        ax3 = axs[2]
        ax4 = axs[3]
        ax5 = axs[4]

        # ax1.set_xlim(-0.5, 3.0)
        # ax2.set_xlim(4.0, 7.0)
        # ax3.set_xlim(8.0, 11.0)
        # ax4.set_xlim(12.0, 15.0)
        # ax5.set_xlim(16.0, 19.6)

        ax1.set_xlim(0, 5)
        ax2.set_xlim(6, 11)
        ax3.set_xlim(12, 17)
        ax4.set_xlim(18, 23)
        ax5.set_xlim(24, 29.6)

        br1 = range(len(df['Total']))
        # br2 = [x + barWidth for x in br1]
        # br3 = [x + barWidth for x in br2]
        # br4 = [x + barWidth for x in br3]

        ax1.errorbar(br1, list(df['COVID-19'].values), yerr=list(df_sem['COVID-19'].values), color='green',
                     linestyle='dashed', linewidth=2, elinewidth=2, label='COVID-19 publications')
        ax2.errorbar(br1, list(df['COVID-19'].values), yerr=list(df_sem['COVID-19'].values), color='green',
                     linestyle='dashed', linewidth=2, elinewidth=2)
        ax3.errorbar(br1, list(df['COVID-19'].values), yerr=list(df_sem['COVID-19'].values), color='green',
                     linestyle='dashed', linewidth=2, elinewidth=2)
        ax4.errorbar(br1, list(df['COVID-19'].values), yerr=list(df_sem['COVID-19'].values), color='green',
                     linestyle='dashed', linewidth=2, elinewidth=2)
        ax5.errorbar(br1, list(df['COVID-19'].values), yerr=list(df_sem['COVID-19'].values), color='green',
                     linestyle='dashed', linewidth=2, elinewidth=2)

        ax1.errorbar(br1, list(df['Non COVID-19'].values), yerr=list(df_sem['Non COVID-19'].values), color='magenta',
                     linestyle='dashed', linewidth=2, elinewidth=2, label='Non COVID-19 publications')
        ax2.errorbar(br1, list(df['Non COVID-19'].values), yerr=list(df_sem['Non COVID-19'].values), color='magenta',
                     linestyle='dashed', linewidth=2, elinewidth=2)
        ax3.errorbar(br1, list(df['Non COVID-19'].values), yerr=list(df_sem['Non COVID-19'].values), color='magenta',
                     linestyle='dashed', linewidth=2, elinewidth=2)
        ax4.errorbar(br1, list(df['Non COVID-19'].values), yerr=list(df_sem['Non COVID-19'].values), color='magenta',
                     linestyle='dashed', linewidth=2, elinewidth=2)
        ax5.errorbar(br1, list(df['Non COVID-19'].values), yerr=list(df_sem['Non COVID-19'].values), color='magenta',
                     linestyle='dashed', linewidth=2, elinewidth=2)

        ax1.errorbar(br1, list(df['Total'].values), yerr=list(df_sem['Total'].values), color='blue', linestyle='dashed', linewidth=2, elinewidth=2,
                     label='Total publications')
        ax2.errorbar(br1, list(df['Total'].values), yerr=list(df_sem['Total'].values), color='blue', linestyle='dashed', linewidth=2, elinewidth=2)
        ax3.errorbar(br1, list(df['Total'].values), yerr=list(df_sem['Total'].values), color='blue', linestyle='dashed', linewidth=2, elinewidth=2)
        ax4.errorbar(br1, list(df['Total'].values), yerr=list(df_sem['Total'].values), color='blue', linestyle='dashed', linewidth=2, elinewidth=2)
        ax5.errorbar(br1, list(df['Total'].values), yerr=list(df_sem['Total'].values), color='blue', linestyle='dashed', linewidth=2, elinewidth=2)


        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        ax3.spines['top'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.spines['left'].set_visible(False)
        ax4.spines['top'].set_visible(False)
        ax5.spines['right'].set_visible(False)
        ax5.spines['left'].set_visible(False)
        ax5.spines['top'].set_visible(False)
        ax1.yaxis.tick_left()

        ax1.tick_params(axis='y', color='white', labelsize='16')
        ax2.tick_params(axis='y', color='white')
        ax3.tick_params(axis='y', color='white')
        ax4.tick_params(axis='y', color='white')
        ax5.tick_params(axis='y', color='white')

        # ax1.set_xticklabels(list(df.T.keys()), rotation=70)
        # ax2.set_xticklabels(list(df.T.keys()[4:8]), rotation=70)
        # ax3.set_xticklabels(list(df.T.keys()[8:12]), rotation=70)
        # ax4.set_xticklabels(list(df.T.keys()[12:16]), rotation=70)
        # ax5.set_xticklabels(list(df.T.keys()[16:]), rotation=70)
        ax1.set_xticklabels(list(df.T.keys()[0:6]), rotation=70, fontdict={'fontsize': 16})
        ax2.set_xticklabels(list(df.T.keys()[6:12]), rotation=70, fontdict={'fontsize': 16})
        ax3.set_xticklabels(list(df.T.keys()[12:18]), rotation=70, fontdict={'fontsize': 16})
        ax4.set_xticklabels(list(df.T.keys()[18:24]), rotation=70, fontdict={'fontsize': 16})
        ax5.set_xticklabels(list(df.T.keys()[24:30]), rotation=70, fontdict={'fontsize': 16})


        d = .005  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (-d, +d), **kwargs)
        ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)

        kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
        ax3.plot((-d, +d), (-d, +d), **kwargs)
        ax3.plot((1 - d, 1 + d), (-d, +d), **kwargs)

        kwargs.update(transform=ax4.transAxes)  # switch to the bottom axes
        ax4.plot((-d, +d), (-d, +d), **kwargs)
        ax4.plot((1 - d, 1 + d), (-d, +d), **kwargs)

        kwargs.update(transform=ax5.transAxes)  # switch to the bottom axes
        ax5.plot((-d, +d), (-d, +d), **kwargs)

        ax3.set_xlabel('First six months in 2016-2020', fontdict={'fontsize': 16})
        ax1.set_ylabel('Average acceptance time (in days)', fontdict={'fontsize': 16})

        fig.legend(loc=(0.07, 0.17), borderaxespad=0.1, fontsize=14)
        ax3.set_title(title, color='black', fontdict={'fontsize': 18})

        plt.show()

    def plt_covid_by_non_covid_acc_time_old(self,df, avgs_for_journals_history, title):
        plt.title(title, color='black')
        month_non_covid_jan = []
        month_covid_jan = []
        labels=[]
        y_pos=0
        # for idx, row in enumerate(df.iterrows()):
        #     print(row)
        #     pub_name = row[0]
        #     total_jan = float(row[1].Total_acc_time_meanApr)
        #     covid_jan = float(row[1].COVID_acc_time_mean_Apr)
        #     total_feb = float(row[1].Total_acc_time_meanFeb)
        #     covid_feb = float(row[1].COVID_acc_time_mean_Feb)
        #     labels.append(pub_name)
        #     month_totals_jan.append(total_jan)
        #     month_covid_jan.append(covid_jan)
            # labels = ['\n'.join(wrap(l, 25)) for l in labels]
            # width = 0.35  # the width of the bars: can also be len(x) sequence
        labels=list(df.index)
        history=list(avgs_for_journals_history)
        month_covid_jan=list(df['COVID_acc_time_mean_Jan'])
        month_non_covid_jan=list(df['non_COVID_acc_time_mean_Jan'])
        month_covid_feb = list(df['COVID_acc_time_mean_Feb'])
        month_non_covid_feb = list(df['non_COVID_acc_time_mean_Feb'])
        month_covid_mar = list(df['COVID_acc_time_mean_Mar'])
        month_non_covid_mar = list(df['non_COVID_acc_time_mean_Mar'])
        month_covid_apr = list(df['COVID_acc_time_mean_Apr'])
        month_non_covid_apr = list(df['non_COVID_acc_time_mean_Apr'])
        y_pos = np.arange(len(labels))
        plt.barh(y_pos+0.15, history, height=0.1, align='center', color='magenta', label='Average 2016-2020')
        plt.barh(y_pos, month_non_covid_jan, height=0.1, align='center', color='green', label='Non COVID')
        plt.barh(y_pos, month_covid_jan, height=0.1, align='center', color='blue', label='COVID')
        plt.barh(y_pos-0.15, month_non_covid_feb, height=0.1, align='center', color='green')
        plt.barh(y_pos-0.15, month_covid_feb, height=0.1, align='center', color='blue' )
        plt.barh(y_pos - 0.3, month_non_covid_mar, height=0.1, align='center', color='green' )
        plt.barh(y_pos - 0.3, month_covid_mar, height=0.1, align='center', color='blue')
        plt.barh(y_pos - 0.45, month_non_covid_apr, height=0.1, align='center', color='green' )
        plt.barh(y_pos - 0.45, month_covid_apr, height=0.1, align='center', color='blue')
            # plt.bar(y_pos+0.2, month_covid, 0.2, align='center', color='blue', label='covid-19 related publications')

        plt.yticks(y_pos, labels, rotation=0, fontsize=8)
        plt.xlabel('Time to Acceptance (in days)')
            # plt.barh(y_pos,month_non_covid, 0.30, color='red', label='total publication')
            # plt.barh(y_pos,month_covid, 0.30, alpha=0.5, color='blue', label='covid-19 related publications')
        plt.legend()
            # plt.yticks(y_pos,labels)
            # plt.xlabel('# of publications')

        plt.show()

    def plt_covid_by_non_covid_acc_time(self, df, avgs_for_journals_history, title):
        plt.title(title, color='black', fontdict={'fontsize': 16})

        labels = list(df.index)
        history = list(avgs_for_journals_history)
        month_covid_jan = list(df['2020-01_COVID'])
        month_non_covid_jan = list(df['2020-01_Non_COVID'])
        month_covid_feb = list(df['2020-02_COVID'])
        month_non_covid_feb = list(df['2020-02_Non_COVID'])
        month_covid_mar = list(df['2020-03_COVID'])
        month_non_covid_mar = list(df['2020-03_Non_COVID'])
        month_covid_apr = list(df['2020-04_COVID'])
        month_non_covid_apr = list(df['2020-04_Non_COVID'])
        month_covid_may = list(df['2020-05_COVID'])
        month_non_covid_may = list(df['2020-05_Non_COVID'])
        month_covid_jun = list(df['2020-06_COVID'])
        month_non_covid_jun = list(df['2020-06_Non_COVID'])
        y_pos = np.arange(len(labels))
        plt.barh(y_pos + 0.3, history, height=0.07, align='center', color='magenta', label='Average 2016-2019')
        plt.barh(y_pos + 0.2, month_non_covid_jan, height=0.07, align='center', color='blue', label='Non COVID')
        plt.barh(y_pos + 0.2, month_covid_jan, height=0.07, align='center', color='green', label='COVID')
        plt.barh(y_pos + 0.1, month_non_covid_feb, height=0.07, align='center', color='blue')
        plt.barh(y_pos + 0.1, month_covid_feb, height=0.07, align='center', color='green')
        plt.barh(y_pos, month_non_covid_mar, height=0.07, align='center', color='blue')
        plt.barh(y_pos, month_covid_mar, height=0.07, align='center', color='green')
        plt.barh(y_pos - 0.1, month_non_covid_apr, height=0.07, align='center', color='blue')
        plt.barh(y_pos - 0.1, month_covid_apr, height=0.07, align='center', color='green')
        plt.barh(y_pos - 0.2, month_non_covid_may, height=0.07, align='center', color='blue')
        plt.barh(y_pos - 0.2, month_covid_may, height=0.07, align='center', color='green')
        plt.barh(y_pos - 0.3, month_non_covid_jun, height=0.07, align='center', color='blue')
        plt.barh(y_pos - 0.3, month_covid_jun, height=0.07, align='center', color='green')
        # plt.bar(y_pos+0.2, month_covid, 0.2, align='center', color='blue', label='covid-19 related publications')

        plt.yticks(y_pos, labels, rotation=0, fontsize=12, wrap=True)
        plt.xlabel('Time to Acceptance (in days)', fontdict={'fontsize': 16})
        # plt.barh(y_pos,month_non_covid, 0.30, color='red', label='total publication')
        # plt.barh(y_pos,month_covid, 0.30, alpha=0.5, color='blue', label='covid-19 related publications')
        plt.legend(loc=(0.7, 0.07),fontsize=14)
        # plt.yticks(y_pos,labels)
        # plt.xlabel('# of publications')

        plt.show()

    def plt_journals_by_month(self,publications, month, month_name):
        labels=[]
        month_totals=[]
        month_covid=[]
        month_in_press=[]
        year_month_str='2020/'+str(month)
        for pub_name, data in publications.items():
            if not year_month_str in data.keys():
                continue
            print(pub_name)
            total = int(data[year_month_str]['total'])
            covid = int(data[year_month_str]['covid'])
            articles_in_press=int(data[year_month_str]['articles_in_press'])
            if int(total)==0:
                continue
            labels.append(pub_name)
            month_totals.append(total)
            month_covid.append(covid)
            month_in_press.append(articles_in_press)
        labels = ['\n'.join(wrap(l, 25)) for l in labels]
        width = 0.35  # the width of the bars: can also be len(x) sequence
        y_pos=np.arange(len(labels))
        plt.barh(y_pos,month_totals, 0.3, align='center', color='red', label='total publications')
        plt.barh(y_pos, month_in_press, 0.3, align='center', color='green', label='articles in press')
        plt.barh(y_pos,month_covid, 0.3, align='center', color='blue', label='covid-19 related publications')
        # plt.bar(y_pos+0.2, month_covid, 0.2, align='center', color='blue', label='covid-19 related publications')

        plt.yticks(y_pos,labels, rotation=0, fontsize=8)
        plt.xlabel('# of publications')
        # plt.barh(y_pos,month_totals, 0.30, color='red', label='total publication')
        # plt.barh(y_pos,month_covid, 0.30, alpha=0.5, color='blue', label='covid-19 related publications')
        plt.legend()
        # plt.yticks(y_pos,labels)
        # plt.xlabel('# of publications')
        plt.title('publication per journal in '+month_name)


        plt.show()

        # plt.savefig('survey.pdf', bbox_inches='tight', pad_inches=0.001, dpi=300)

    def plt_country_by_month(self, month_dict,title ):
        x=month_dict.keys()
        y=month_dict.values()
        plt.figure()
        plt.plot(list(x)[0:15],list(y)[0:15])
        plt.xticks(rotation=30, ha='right')
        plt.title(title)
        plt.show()

    def plt_countries_totals_history(self, year_df):
        xs,ys=year_df[1].values, year_df[0].values
        ys = list(np.around(np.array(ys), 2))
        plt.scatter(xs,ys)
        idx=0
        for x, y in zip(xs, ys):
            label = "{} ({})".format(year_df[0].index[idx], np.round(year_df[0].values[idx],2))
            xy_coord=(-4, 2)
            if x==2016:
                xy_coord=(1,1)
                if idx==0:
                    xy_coord = (0, 10)
                if idx==1 or idx==4:
                    xy_coord = (0, -13)
                plt.annotate(label,  # this is the text
                             (x, y),  # this is the point to label
                             textcoords="offset points",  # how to position the text
                             xytext=xy_coord,  # distance from text to points (x,y)
                             ha='left',
                             fontsize=14)  # horizontal alignment can be left, right or center
            else:
                if (x == 2017 or x==2018 or x==2019) and idx==4:
                        xy_coord = (-4, -7)
                plt.annotate(label,  # this is the text
                         (x, y),  # this is the point to label
                         textcoords="offset points",  # how to position the text
                         xytext=xy_coord,  # distance from text to points (x,y)
                         ha='right',
                         fontsize=14 )  # horizontal alignment can be left, right or center
            idx+=1
        # plt.xticks([2016,2017,2018,2019,2020,2022])
        plt.gca().set_xticklabels(['','2016', '2017', '2018', '2019', '2020-cov', '2020-non-cov'], fontdict={'fontsize': 16})
        plt.ylabel('Average number of articles (over six months)', fontdict={'fontsize': 16})
        plt.tick_params(axis='y',labelsize='16')
        # plt.show()


    def plt_countries_collab_diversity_history(self, year_df):

        xs,ys=year_df['year'].values, year_df['num_collaborating_countries'].values
        plt.scatter(xs,ys)
        idx=0
        for x, y in zip(xs, ys):
            label = "{} ({})".format(year_df['country'].values[idx],year_df['num_collaborating_countries'].values[idx])
            xy_coord=(-3, 2)
            if x==2016:
                xy_coord=(3,1)
                if idx==3:
                    xy_coord = (3, -6)
                if idx==4:
                    xy_coord = (3, -20)
                plt.annotate(label,  # this is the text
                             (x, y),  # this is the point to label
                             textcoords="offset points",  # how to position the text
                             xytext=xy_coord,  # distance from text to points (x,y)
                             ha='left',
                             fontsize=14)  # horizontal alignment can be left, right or center
            else:
                if (x==2018 or x==2019 or x==2020) and idx==4:
                        xy_coord = (-3, -12)
                if x==2019 and idx==1:
                    xy_coord = (-3, -8)
                # if x==2020 and idx==3:
                #     xy_coord = (-5, -2)
                if x==2020 and idx==1:
                    xy_coord = (-3, -6)
                plt.annotate(label,  # this is the text
                         (x, y),  # this is the point to label
                         textcoords="offset points",  # how to position the text
                         xytext=xy_coord,  # distance from text to points (x,y)
                         ha='right',
                         fontsize=14 )  # horizontal alignment can be left, right or center
            idx+=1
        plt.xticks([2016,2017,2018,2019,2020,2021], fontsize=16)
        plt.gca().set_xticklabels(['2016', '2017', '2018', '2019', '2020-cov', '2020-non-cov'], fontdict={'fontsize': 16})
        plt.ylabel('Number of collaborating countries', fontdict={'fontsize': 16})
        # plt.show()

    def plt_covid_countries_collab_diversity(self, month_df, color):
        xs, ys = month_df['Month'].values, month_df['num_collaborating_countries'].values
        plt.scatter(xs, ys, color=color)
        month=xs[0]
        idx = 0
        prev_y=500
        prev_y_coord=50
        x_coord = -5
        y_coord = 2
        ha='right'
        for x, y in zip(xs, ys):
            label = "{} ({})".format(month_df['country'].values[idx], month_df['num_collaborating_countries'].values[idx])
            if x==1:
                x_coord=5
                ha='left'
            if prev_y-y<=4:
                y_coord=prev_y_coord-4
            if prev_y==y:
                y_coord=prev_y_coord-9
            if x==5 and color=='green':
                    if (idx==2 ):
                        y_coord=-12
                    if idx==3:
                        y_coord=-18
                    if idx==4:
                        y_coord=-26



            xy_coord=(x_coord,y_coord)
            plt.annotate(label,  # this is the text
                         (x, y),  # this is the point to label
                         textcoords="offset points",  # how to position the text
                         xytext=xy_coord,  # distance from text to points (x,y)
                         ha=ha,
                         fontsize=12)  # horizontal alignment can be left, right or center
            idx+=1
            prev_y=y
            prev_y_coord=y_coord
        plt.gca().set_xticklabels(['','2020-01', '2020-02', '2020-03', '2020-04', '2020-05', '2020-06'], rotation=45,  fontdict={'fontsize': 16})
        plt.xlabel('First six months in 2020', fontdict={'fontsize': 16})
        plt.ylabel('Number of collaborating countries', fontdict={'fontsize': 16})


    def plt_countries_collab_num_papers_history(self, year_df):

        xs,ys=year_df['year'].values, year_df['avg_num_papers'].values
        plt.scatter(xs,ys)
        idx=0
        prev_y = 500
        prev_y_coord = 50
        x_coord = 0
        y_coord = 2
        ha = 'right'
        for x, y in zip(xs, ys):
            label = "{} ({}%, {})".format(year_df['country'].values[idx],np.round(year_df['avg_num_papers'].values[idx],2),year_df['num_collaborate_papers'].values[idx])
            if x==2016:
                x_coord = -55
                ha = 'left'
                if idx==1:
                    y_coord=-16
                if idx==2:
                    y_coord=16
                if idx==3:
                    y_coord=-2
                if idx==4:
                    y_coord=-25
                #     xy_coord = (5, -4)
                # if idx==4:
                #     xy_coord = (5, -10)

                # if (x==2018 or x==2019 or x==2020) and idx==4:
                #         xy_coord = (-5, -4)
                # if x==2019 and idx==1:
                #     xy_coord = (-5, -2)
                # # if x==2020 and idx==3:
                # #     xy_coord = (-5, -2)
                # if x==2020 and idx==1:
                #     xy_coord = (-5, -2)
            else:
                if x == 2017:
                    if idx == 2 or idx == 3 or idx==4:
                        y_coord = 18
                if prev_y-y<=4:
                    y_coord=prev_y_coord-6
                if prev_y==y:
                    y_coord=prev_y_coord-10

            xy_coord=(x_coord,y_coord)
            plt.annotate(label,  # this is the text
                    (x, y),  # this is the point to label
                    textcoords="offset points",  # how to position the text
                    xytext=xy_coord,  # distance from text to points (x,y)
                    ha=ha,
                    fontsize=14 )  # horizontal alignment can be left, right or center
            prev_y = y
            prev_y_coord = y_coord
            idx+=1
        plt.xticks([2016,2017,2018,2019,2020,2021], fontsize=16)
        plt.gca().set_xticklabels(['2016', '2017', '2018', '2019', '2020-cov', '2020-non-cov'], fontdict={'fontsize': 16})
        plt.ylabel('collaboration papers percentage', fontdict={'fontsize': 16})
        # plt.show()

    def plt_covid_countries_collab_num_papers(self, month_df, color, perc=False):
        xs, ys = month_df['Month'].values, month_df['num_collaborate_papers'].values
        if perc:
            xs, ys = month_df['Month'].values, month_df['avg_num_papers'].values
        plt.scatter(xs, ys, color=color)
        idx = 0
        prev_y=500
        prev_y_coord=50
        x_coord = -3
        y_coord = 2
        ha='right'
        for x, y in zip(xs, ys):
            label = "{} ({})".format(month_df['country'].values[idx],month_df['num_collaborate_papers'].values[idx])
            if x==1:
                x_coord=1
                y_coord=14
                ha='left'
                if idx==3:
                    x_coord=-24
            # if x==3 and idx==2 and color=='green':
            #     plt.scatter(x,y+4, color=color)
            #     y_coord=20
            # if x==5 and idx==2 and color=='blue':
            #     y_coord=-2
            # if x==5 and idx==2 and color=='green':
            #     plt.scatter(x,y+2, color=color)
            #     y_coord=10
            if x==2 and idx==3:
                y_coord=-12
            if x==3 and color=='green':
                y_coord=-6
            if x==5 and color=='green' and (idx==2 or idx==3 or idx==4):
                y_coord=-20
            if x == 5 and color == 'blue' and (idx==2 or idx == 3 or idx == 4):
                y_coord = 20
            if prev_y-y<=4:
                y_coord=prev_y_coord-6
            if prev_y==y:
                y_coord=prev_y_coord-12


            # if x==5 and idx==4 and color=='green':
            #     y_coord=-10
            # if x == 3 and idx == 0 and color == 'green':
            #     plt.scatter(x,y+1, color=color)
            #     y_coord = 10
            # if x == 3 and idx == 1 and color == 'green':
            #     y_coord = -2

            xy_coord=(x_coord,y_coord)
            plt.annotate(label,  # this is the text
                         (x, y),  # this is the point to label
                         textcoords="offset points",  # how to position the text
                         xytext=xy_coord,  # distance from text to points (x,y)
                         ha=ha,
                         fontsize=12)  # horizontal alignment can be left, right or center
            idx+=1
            prev_y=y
            prev_y_coord=y_coord
        plt.gca().set_xticklabels(['','2020-01', '2020-02', '2020-03', '2020-04', '2020-05', '2020-06'], rotation=45, fontdict={'fontsize': 16})
        plt.xlabel('First six months in 2020', fontdict={'fontsize': 16})
        if perc:
            plt.ylabel('Collaboration papers percentage', fontdict={'fontsize': 16})
        else:
            plt.ylabel('Number of collaboration papers', fontdict={'fontsize': 16})


    def plt_scatter(self,df):
        df.plot.scatter(x='SJR', y='sum_COVID')




    def show(self, title, legend=None, color=None, fontdict=None):
        if legend:
            plt.legend(legend, fontsize=14)
        plt.title(title, color=color, fontdict=fontdict)
        plt.show()

    def plt_num_authors_by_month(self,month_name,med,bio, arxiv):
        y_pos=np.arange(15)
        plt.bar(y_pos,med*100, 0.2, align='center', color='blue', label='medrxiv')
        plt.bar(y_pos+0.2, bio*100, 0.2, align='center', color='magenta', label='biorxiv')
        plt.bar(y_pos+0.4, arxiv*100, 0.2, align='center', color='green', label='arxiv')
        plt.ylabel('% of articles')
        plt.xlabel('num authors, 14 includes 14 or more authors')
        plt.legend()
        plt.title('num authors in article in '+month_name)
        plt.show()

    def plt_acceptance_speed(self,data_dict,title):
        x = data_dict.keys()
        vals = data_dict.values()
        totals=[]
        week=[]
        fortnight=[]
        triweek=[]
        month_plus=[]
        for k, item in data_dict.items():
            totals.append(item['Total_days'])
            week.append(item['acc_speed']['week'])
            fortnight.append(item['acc_speed']['fortnight'])
            triweek.append(item['acc_speed']['triweek'])
            month_plus.append(item['acc_speed']['month_plus'])

        plt.figure()
        plt.plot(list(x), totals, color='blue', linestyle="-", label='totals')
        plt.plot(list(x), week, color='green', linestyle=":", label='week')
        plt.plot(list(x), fortnight, color='magenta', linestyle="--", label='fortnight')
        plt.plot(list(x), triweek, color='red', linestyle="--", label='triweek')
        plt.plot(list(x), month_plus, color='cyan', linestyle="--", label='month plus')

        plt.legend()
        plt.xticks(rotation=30, ha='right')
        plt.title(title)
        plt.show()


if __name__ == '__main__':
    vis=Visualization()
    pss=ParseScopusSearch()
    scopus_path = 'D:\\shir\\study\\covid_19\\scopus'
    utils=Utils(path=scopus_path)
    obj=utils.load_obj('journals_data_apr')
    name='journals_data_inc_apr.csv'
    file_path = os.path.join(scopus_path, name)
    with open(file_path, 'w') as f:  # Just use 'w' mode in 3.x
        w = csv.DictWriter(f, obj.keys())
        w.writeheader()
        w.writerow(obj)
    vis.plt_journals_by_month(obj,1,"January")
    vis.plt_journals_by_month(obj,2,"February")
    vis.plt_journals_by_month(obj,3,"March")
    vis.plt_journals_by_month(obj,4,"April")


    medrxiv_data=[[0,0,0,2,120,533,1265],[162,201,170,226,433,805,1667],[0,0,0,0.008849558,0.277136259,0.662111801,0.75884823]]
    bioarxiv_data=[[0,0,0,12,54,138,348],[3111,2801,2636,2945,3248,3584,4230],[0,0,0, 0.004074703,0.016625616, 0.038504464,0.082269504]]
    arxiv_data=[[0,0,0,1,24,126,245],[265,298,240,206, 272,358,492],[0,0,0,0, 0.088235294, 0.351955307, 0.49796748]]
    scopus_data=[[0,0,0,18,182,680,2866]]
    world_infectious_statistics=[[0,0.001,0.266,9.826,85.403,750.890,3090.445]]

    # vis.plt_by_month(medrxiv_data,bioarxiv_data,arxiv_data, scopus_data, world_infectious_statistics)