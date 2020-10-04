import json
from utils import Utils
import os
from AuthorsAndCountries import AuthorsAndCountries
import numpy as np
from visualization import Visualization
import arxiv
import pandas as pd



path='D:\\shir\\study\\covid_19\\'
max_num_authors=15

class ParseArxivs:

    def __init__(self):

        self.articles_by_num_authors_per_month = dict()
        self.articles_by_num_authors_per_month['med']=dict()
        self.articles_by_num_authors_per_month['med']['Jan']=np.zeros(max_num_authors)
        self.articles_by_num_authors_per_month['med']['Feb']=np.zeros(max_num_authors)
        self.articles_by_num_authors_per_month['med']['Mar']=np.zeros(max_num_authors)
        self.articles_by_num_authors_per_month['med']['Apr']=np.zeros(max_num_authors)

        self.articles_by_num_authors_per_month['bio']=dict()
        self.articles_by_num_authors_per_month['bio']['Jan'] = np.zeros(max_num_authors)
        self.articles_by_num_authors_per_month['bio']['Feb'] = np.zeros(max_num_authors)
        self.articles_by_num_authors_per_month['bio']['Mar'] = np.zeros(max_num_authors)
        self.articles_by_num_authors_per_month['bio']['Apr'] = np.zeros(max_num_authors)

        self.articles_by_num_authors_per_month['arxiv'] = dict()
        self.articles_by_num_authors_per_month['arxiv']['Jan'] = np.zeros(max_num_authors)
        self.articles_by_num_authors_per_month['arxiv']['Feb'] = np.zeros(max_num_authors)
        self.articles_by_num_authors_per_month['arxiv']['Mar'] = np.zeros(max_num_authors)
        self.articles_by_num_authors_per_month['arxiv']['Apr'] = np.zeros(max_num_authors)


    def load_data(self, file):
        with open(file, 'r', encoding='utf8') as myfile:
            data = myfile.read()
        response_json = json.loads(data)
        # print(response_json)
        return response_json

    def load_arxiv(self):
        query_vals="(ti:covid OR ti:nCoV OR ti:coronavirus OR abs:covid  OR abs:nCov OR abs:coronavirus)"
        query_categories="(cat:q-bio.BM OR cat:q-bio.CB OR cat:q-bio.GN OR cat:q-bio.MN OR cat:q-bio.NC OR cat:q-bio.OT OR cat:q-bio.PE OR cat:q-bio.QM OR cat:q-bio.SC OR cat:q-bio.TO)"
        result = arxiv.query(query=query_vals, max_results=1000, iterative=False)
        idx=0
        filtered_result=[]
        for paper in result:
            if paper['published_parsed'][1]==5 or paper['published_parsed'][0]!=2020:
                continue
            if paper['title'].lower().find('covid-19')==-1 and  paper['title'].lower().find('coronavirus')==-1 and paper['title'].lower().find('ncov')==-1:
                if paper['summary'].lower().find('covid-19') == -1 and paper['summary'].lower().find('coronavirus') == -1 and paper['summary'].lower().find('ncov') == -1:
                    continue
            # if paper['arxiv_primary_category']['term'].lower().find('q-bio') == -1:
            #     foundBio=False
            #     for tag in paper['tags']:
            #         if tag['term'].lower().find('q-bio') != -1:
            #             foundBio=True
            #         if foundBio:
            #             continue
            #     if not foundBio:
            #         print("cat is {}".format(paper['arxiv_primary_category']['term']))
            #         print(paper)
            #         continue
            idx+=1
            filtered_result.append(paper)
        print(len(result))
        print(idx)
        return filtered_result

    def split_arxiv_by_month(self, results):
        arxiv = dict()

        arxiv['Jan'] = []
        arxiv['Feb'] = []
        arxiv['Mar'] = []
        arxiv['Apr'] = []
        for paper in results:
            if paper['published_parsed'][1]==1:
                arxiv['Jan'].append(paper)
            if paper['published_parsed'][1] == 2:
                arxiv['Feb'].append(paper)
            if paper['published_parsed'][1] == 3:
                arxiv['Mar'].append(paper)
            if paper['published_parsed'][1] == 4:
                arxiv['Apr'].append(paper)
        return arxiv

    def split_by_month(self,json):
        print(type(json))
        medrvix=dict()
        biorxiv=dict()

        medrvix['Jan']=[]
        medrvix['Feb']=[]
        medrvix['Mar']=[]
        medrvix['Apr']=[]

        biorxiv['Jan'] = []
        biorxiv['Feb'] = []
        biorxiv['Mar'] = []
        biorxiv['Apr'] = []

        for item in json['rels']:
            if str(item['rel_date']).startswith('2020-04'):
                if str(item['rel_site']).startswith('medrxiv'):
                    medrvix['Apr'].append(item)
                else:
                    biorxiv['Apr'].append(item)
            if str(item['rel_date']).startswith('2020-03'):
                if str(item['rel_site']).startswith('medrxiv'):
                    medrvix['Mar'].append(item)
                else:
                    biorxiv['Mar'].append(item)
            if str(item['rel_date']).startswith('2020-02'):
                if str(item['rel_site']).startswith('medrxiv'):
                    medrvix['Feb'].append(item)
                else:
                    biorxiv['Feb'].append(item)
            if str(item['rel_date']).startswith('2020-01'):
                if str(item['rel_site']).startswith('medrxiv'):
                    medrvix['Jan'].append(item)
                else:
                    biorxiv['Jan'].append(item)


        return medrvix,biorxiv

    # def extract_by_doi(self,doi):

    def get_longtitudal_num_papers_changes(self):
        data = utils.load_csv_data_to_df('covid_counts_preprints.csv')
        data.set_index('Venue', inplace=True)
        data = data.apply(pd.to_numeric, errors='coerce')
        data=data.T
        data['world_stats']=data['world_stats']/5000
        vis.plt_rxivs_cov_and_totals(data, "COVID-19 and Total publications, Pre-print servers")



if __name__ == '__main__':
    pa=ParseArxivs()
    utils=Utils(path=path)
    vis=Visualization()



    pa.get_longtitudal_num_papers_changes()
    exit(0)

    results=pa.load_arxiv()
    arxiv=pa.split_arxiv_by_month(results)
    print('Apr arxiv {}'.format(len(arxiv['Apr'])))
    print('Mar arxiv {}'.format(len(arxiv['Mar'])))
    print('Feb arxiv {}'.format(len(arxiv['Feb'])))
    print('Jan arxiv {}'.format(len(arxiv['Jan'])))

    print(arxiv['Jan'])

    res= pa.load_data(os.path.join(path,'collection_json.json'))
    med,bio=pa.split_by_month(res)
    print('Apr med {}, bio {}'.format(len(med['Apr']),len(bio['Apr'])))
    print('Mar med {}, bio {}'.format(len(med['Mar']),len(bio['Mar'])))
    print('Feb med {}, bio {}'.format(len(med['Feb']),len(bio['Feb'])))
    print('Jan med {}, bio {}'.format(len(med['Jan']),len(bio['Jan'])))
    anc=AuthorsAndCountries()
    vis = Visualization()
    for key in med.keys():
        curr_month=key
        for article in med[curr_month]:
            doi=article['rel_doi']
            num_authors=article['rel_num_authors']
            # res=anc.get_authors_data_by_doi(doi)
            # if isinstance(res, list):
            #     print(res)
            if num_authors>max_num_authors-1:
                # print(num_authors)
                num_authors=max_num_authors-1
            pa.articles_by_num_authors_per_month['med'][curr_month][num_authors]+=1
        for article in bio[curr_month]:
            doi=article['rel_doi']
            num_authors=article['rel_num_authors']
            # res=anc.get_authors_data_by_doi(doi)
            # if isinstance(res, list):
            #     print(res)
            if num_authors>max_num_authors-1:
                # print(num_authors)
                num_authors=max_num_authors-1
            pa.articles_by_num_authors_per_month['bio'][curr_month][num_authors]+=1
        for article in arxiv[curr_month]:
            doi=article['doi']
            if doi!=None:
                res=anc.get_authors_data_by_doi(doi)
                if isinstance(res, list):
                    print(res)
            if article['journal_reference']!=None:
                print("in journal")
            num_authors = len(article['authors'])
            if num_authors>max_num_authors-1:
                # print(num_authors)
                num_authors=max_num_authors-1
            pa.articles_by_num_authors_per_month['arxiv'][curr_month][num_authors]+=1
        print(pa.articles_by_num_authors_per_month['arxiv'][curr_month])
        vis.plt_num_authors_by_month(curr_month,pa.articles_by_num_authors_per_month['med'][curr_month]/len(med[curr_month]),pa.articles_by_num_authors_per_month['bio'][curr_month]/len(bio[curr_month]), pa.articles_by_num_authors_per_month['arxiv'][curr_month]/len(arxiv[curr_month]))