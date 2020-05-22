
import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap
from utils import Utils
import math
from tqdm import tqdm_notebook
# %matplotlib inline
from parse_scopus_search import ParseScopusSearch

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
        labels = ['\n'.join(wrap(l, 10)) for l in labels]
        width = 0.35  # the width of the bars: can also be len(x) sequence
        y_pos=np.arange(len(labels))
        plt.bar(y_pos,month_totals, 0.2, align='center', color='red', label='total publications')
        plt.bar(y_pos, month_in_press, 0.2, align='center', color='green', label='articles in press')
        plt.bar(y_pos,month_covid, 0.2, align='center', color='blue', label='covid-19 related publications')
        # plt.bar(y_pos+0.2, month_covid, 0.2, align='center', color='blue', label='covid-19 related publications')

        plt.xticks(y_pos,labels, rotation=90)
        plt.ylabel('# of publications')
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
            totals.append(item['Total'])
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
    utils=Utils(path='D:\\shir\\study\\covid_19\\scopus')
    obj=utils.load_obj('journals_data_apr')
    vis.plt_journals_by_month(obj,1,"January")
    vis.plt_journals_by_month(obj,2,"February")
    vis.plt_journals_by_month(obj,3,"March")
    vis.plt_journals_by_month(obj,4,"April")


    medrxiv_data=[[0,0,0,2,120,533,1265],[162,201,170,226,433,805,1667],[0,0,0,0.008849558,0.277136259,0.662111801,0.75884823]]
    bioarxiv_data=[[0,0,0,12,54,138,348],[3111,2801,2636,2945,3248,3584,4230],[0,0,0, 0.004074703,0.016625616, 0.038504464,0.082269504]]
    arxiv_data=[[0,0,0,1,24,126,245],[295,347,300,206, 272,358,492],[0,0,0,0, 0.069970845, 0.318987342, 0.690031153]]
    scopus_data=[[0,0,0,19,187,684,2906]]
    world_infectious_statistics=[[0,0.001,0.266,9.826,85.403,750.890,3090.445]]

    vis.plt_by_month(medrxiv_data,bioarxiv_data,arxiv_data, scopus_data, world_infectious_statistics)