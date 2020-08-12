import requests
import json
import pandas as pd
import os
import time
from bio import Entrez
import pickle

path='D:\\shir\\study\\covid_19\\scopus'
class ParseJournalsMetrics:

    headers = {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'x-els-apikey': 'fee04e1e02e2b7ed9d4991b87f29896f',
    }

    sciencedirect = 'https://api.elsevier.com/content/serial/title'


    def send_request(self,url,headers,data):
        query_str=url+'?title=\"'+data+'\"&field=SJR,SNIP,citeScoreYearInfoList&view=STANDARD'
        response = requests.get(query_str, headers=headers)
        response_dict = json.loads(response.text)
        print(response)
        return response_dict

    def iterate_journals(self):
        journals = self.read_from_csv('Journals_list.csv').split('\n')
        journals = journals[0:25]
        # journals=['Journal of Vascular Surgery','Tropical Medicine and International Health']
        dict_ = {'journal_name': [], 'SJR': [], 'SJR_year': [], 'SNIP': [], 'SNIP_year':[], 'CiteScore': [], 'CiteScore_year':[], 'url':[]}

        for idx, journal_name in enumerate(journals):
            # journal_name,value=journal.split(',')
            # if idx % (len(journals) // 10) == 0:
            #     print(f'Processing index: {idx} of {len(journals)}')
            dict_['journal_name'].append(journal_name)
            # dict_['count'].append(value)
            response_dict = self.send_request(pj.sciencedirect, pj.headers, data=journal_name)
            if 'serial-metadata-response' in response_dict.keys() and 'entry' in response_dict['serial-metadata-response'].keys():
                metrics=response_dict['serial-metadata-response']['entry'][0]
                if 'SNIPList' in metrics.keys():
                    Snip_val=metrics['SNIPList']['SNIP'][0]['$']
                    Snip_year=metrics['SNIPList']['SNIP'][0]['@year']
                else:
                    Snip_val='unknwon'
                    Snip_year='unknwon'
                if 'SJRList' in metrics.keys():
                    Sjr_val=metrics['SJRList']['SJR'][0]['$']
                    Sjr_year=metrics['SJRList']['SJR'][0]['@year']
                else:
                    Sjr_val='unknwon'
                    Sjr_year='unknwon'
                if 'citeScoreYearInfoList' in metrics.keys():
                    citeScore_val=metrics['citeScoreYearInfoList']['citeScoreCurrentMetric']
                    citeSCore_year=metrics['citeScoreYearInfoList']['citeScoreCurrentMetricYear']
                else:
                    citeScore_val = 'unknwon'
                    citeSCore_year = 'unknown'
                if 'prism:url' in metrics.keys():
                    issn_url=metrics['prism:url']
                else:
                    issn_url="unknown"
            else:
                Snip_val = 'unknwon'
                Snip_year = 'unknwon'
                Sjr_val = 'unknwon'
                Sjr_year = 'unknwon'
                citeScore_val = 'unknwon'
                citeSCore_year = 'unknown'
                issn_url = "unknown"
            dict_['SNIP'].append(Snip_val)
            dict_['SNIP_year'].append(Snip_year)
            dict_['SJR'].append(Sjr_val)
            dict_['SJR_year'].append(Sjr_year)
            dict_['CiteScore'].append(citeScore_val)
            dict_['CiteScore_year'].append(citeSCore_year)
            dict_['url'].append(issn_url)

        df_covid = pd.DataFrame(dict_,
                    columns=['journal_name', 'SJR', 'SJR_year', 'SNIP', 'SNIP_year', 'CiteScore', 'CiteScore_year', 'url'])
        return df_covid



    def write_to_csv(self, df, name):
        df.to_csv(name, index=False)

    def read_from_csv(self, name):
        csv_path = os.path.join(path, name)
        with open(csv_path, encoding="utf8") as csv_file:
            contents = csv_file.read()
        return contents

if __name__ == '__main__':
    pj=ParseJournalsMetrics()
    df=pj.iterate_journals()
    pj.write_to_csv(df,os.path.join(path,'journals_list_metrics_tmp.csv'))
    print(df)