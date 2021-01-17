import requests
import json
import pandas as pd
import os
import time
from Bio import Entrez
from utils import Utils
import pickle

retmax=400
path='D:\\shir\\study\\covid_19\\scopus'
class ParseScopusSearch:

    headers = {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'x-els-apikey': 'fee04e1e02e2b7ed9d4991b87f29896f',
    }

    data = '{\n  "qs": "\\"COVID-19\\" OR \\"ncov\\" AND 2020",\n "display": {\n      "offset": 100,\n      "show": 100,\n      "sortBy": "date"\n  }\n}'
    sciencedirect = 'https://api.elsevier.com/content/search/sciencedirect'


    def send_request(self,url,headers,data):
        response = requests.put(url, headers=headers, data=data)
        response_dict = json.loads(response.text)
        print(response)
        if response.status_code==200:
            return response_dict
        else:
            return None


    def get_breaks(self, content, length):
        data = ""
        words = content.split(' ')
        total_chars = 0

        # add break every length characters
        for i in range(len(words)):
            total_chars += len(words[i])
            if total_chars > length:
                data = data + "<br>" + words[i]
                total_chars = 0
            else:
                data = data + " " + words[i]
        return data


    def convert_to_df(self, all_json):
        dict_ = {'title': [], 'authors': [], 'doi': [], 'journal': [], 'online_date':[], 'uri':[]}

        for index, part in enumerate(all_json):
            for idx, entry in enumerate(part):
                if idx % (len(part) // 10) == 0:
                    print(f'Processing index: {idx} of {len(all_json)}')
                # content = FileReader(entry)

                try:
                    # if more than one author
                    authors = entry['authors']
                    auth=""
                    if len(authors) > 1:
                        for author in authors:
                            auth+=author['name']+";"
                            # more than 2 authors, may be problem when plotting, so take first 2 append with ...
                        dict_['authors'].append(auth)
                    else:
                        # authors will fit in plot
                        dict_['authors'].append(authors[0]['name'])
                except Exception as e:
                    # if only one author - or Null valie
                    dict_['authors'].append("unknown")

                # add the title information, add breaks when needed
                try:
                    title = self.get_breaks(entry['title'], 40)
                    dict_['title'].append(title)

                # if title was not provided
                except Exception as e:
                    dict_['title'].append(entry['title'])

                # add the journal information
                dict_['journal'].append(entry['sourceTitle'])
                dict_['doi'].append(entry['doi'])
                date=entry['loadDate']
                date=date[0:date.index('T')]
                dict_['online_date'].append(date)
                dict_['uri'].append(entry['uri'])


        df_covid = pd.DataFrame(dict_, columns=[ 'title', 'authors','doi', 'journal','online_date', 'uri'])
        df_covid.drop_duplicates(['title'], inplace=True)
        return df_covid

    def get_journal_data(self,journal_name, min_date, max_date):
        query=journal_name+'[journal]'
        # query='Travel+Medicine+and+Infectious+Disease[journal]'
        Entrez.email='shirAviv@protonmail.com'
        handle = Entrez.esearch(db='pubmed',
                                sort='date',
                                retmax=retmax,
                                retmode='xml',
                                mindate=min_date,
                                maxdate=max_date,
                                term=query)
        results = Entrez.read(handle)
        if int(results['Count']) > retmax:
            print('more results for {} exist'.format(journal_name))
        return results


        # results = []
        # offset = 0
        # str_offset = str(offset)
        # data = '{\n  "qs": "\\"'+journal_name+'\\" AND 2020",\n "display": {\n      "offset": '+str_offset+',\n      "show": 100,\n      "sortBy": "date"\n  }\n}'
        # res=self.send_request(self.sciencedirect,self.headers,data)
        # results.append(res['results'])
        # numResults=res['resultsFound']
        # offset+=100
        # while offset<numResults:
        #     time.sleep(1)
        #     str_offset = str(offset)
        #     data = '{\n  "qs": "\\"' + journal_name + '\\" AND 2020",\n "display": {\n      "offset": ' + str_offset + ',\n      "show": 100,\n      "sortBy": "date"\n  }\n}'
        #     res = self.send_request(self.sciencedirect, self.headers, data)
        #     results.append(res['results'])
        #     offset += 100
        return results

    def extract_journal_data(self, ids, journal_name):
        ids = ','.join(ids['IdList'])
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=ids)
        results = Entrez.read(handle)
        return results
        # dict_ = {'journal': [], 'online_date': [], 'uri': []}
        # dates_counter=[0,0,0,0,0,0]
        # for index, part in enumerate(all_json):
        #     for idx, entry in enumerate(part):
        #         if idx % (len(part) // 10) == 0:
        #             print(f'Processing index: {idx} of {len(all_json)}')
        #         # content = FileReader(entry)
        #
        #         try:
        #             # add the journal information
        #             journal=entry['sourceTitle']
        #             if (not journal.lower()==journal_name.lower()):
        #                 continue
        #
        #             date = entry['loadDate']
        #             date = date[0:date.index('T')]
        #             if not (date.startswith('2020') or date.startswith('2019')):
        #                 continue
        #             if date[0:7]=='2020-03':
        #                 dates_counter[0]+=1
        #             if date[0:7] == '2020-02':
        #                 dates_counter[1] += 1
        #             if date[0:7] == '2020-01':
        #                 dates_counter[2] += 1
        #             if date[0:7] == '2019-12':
        #                 dates_counter[3] += 1
        #             if date[0:7] == '2019-11':
        #                 dates_counter[4] += 1
        #             if date[0:7] == '2019-10':
        #                 dates_counter[5] += 1
        #
        #             dict_['online_date'].append(date)
        #             dict_['journal'].append(entry['sourceTitle'])
        #             dict_['uri'].append(entry['uri'])
        #         except Exception as e:
        #             continue
        #
        #
        # print(dates_counter)

    def iterate_journals(self):
        num_publication = dict()
        empty_results=set(())
        journals=utils.read_from_csv('Journals_by_month_keyword_short_2.csv').split('\n')
        journals=journals[0:-1]
        for idx,journal in enumerate(journals):
            if idx==100:
                break
            journal_name,value,articles_in_press=journal.split(',')
            if journal_name.lower().endswith('jan'):
                year = 2020
                month = 1
                continue
            if journal_name.lower().endswith('feb'):
                year = 2020
                month = 2
                continue
            if journal_name.lower().endswith('mar'):
                year = 2020
                month = 3
                continue
            if journal_name.lower().endswith('apr'):
                year = 2020
                month = 4
                continue
            if not journal_name in num_publication.keys():
                num_publication[journal_name]=dict()
            str_year_month = str(year) + '/' + str(month)
            num_publication[journal_name][str_year_month]=dict()
            num_publication[journal_name][str_year_month]['covid'] = value
            num_publication[journal_name][str_year_month]['articles_in_press'] = articles_in_press

            if journal_name.lower().endswith('the lancet'):
                empty_results.add(journal_name)
                continue
            res_journals = pss.get_journal_data(journal_name, str_year_month, str_year_month)
            num_publication[journal_name][str_year_month]['total'] = str(res_journals['Count'])
            if int(res_journals['Count'])==0:
                tmp_journal_name=journal_name.replace("&","")
                res_journals = pss.get_journal_data(tmp_journal_name, str_year_month, str_year_month)
                num_publication[journal_name][str_year_month]['total'] = str(res_journals['Count'])
                if int(res_journals['Count']) == 0:
                    empty_results.add(journal_name+" "+str(month))

            print('working on {}'.format(journal_name))

            year_older = 2019
            month_older = 12
            str_year_month = str(year_older) + '/' + str(month_older)
            if  (not str_year_month in num_publication[journal_name].keys()):
                for i in range(3):
                    num_publication[journal_name][str_year_month] = dict()

                    res_journals = pss.get_journal_data(journal_name, str_year_month, str_year_month)
                    num_publication[journal_name][str_year_month]['total'] = str(res_journals['Count'])
                    month_older -= 1
                    str_year_month = str(year_older) + '/' + str(month_older)

        for k,v in num_publication.items():
            print(k,v)
        print(empty_results)
        print(len(empty_results))
        return num_publication
        # pss.extract_journal_data(res_journals, journal_name)

    def handle_special_journals(self, num_publications):
        journal_name='New Scientist'
        # str_year_month = '2020/1'
        # num_publications[journal_name][str_year_month]['total'] = '184'
        str_year_month = '2020/2'
        num_publications[journal_name][str_year_month]['total'] = '237'
        str_year_month = '2020/3'
        num_publications[journal_name][str_year_month]['total'] = '185'

        journal_name='The Lancet'
        str_year_month = '2020/1'
        num_publications[journal_name][str_year_month]['total'] = '142'
        str_year_month = '2020/2'
        num_publications[journal_name][str_year_month]['total'] = '148'
        str_year_month = '2020/3'
        num_publications[journal_name][str_year_month]['total'] = '154'
        str_year_month = '2020/4'
        num_publications[journal_name][str_year_month]['total'] = '172'

        journal_name='Journal of Microbiology Immunology and Infection'
        # str_year_month = '2020/1'
        # num_publications[journal_name][str_year_month]['total'] = '24'
        # str_year_month = '2020/2'
        # num_publications[journal_name][str_year_month]['total'] = '43'
        str_year_month = '2020/3'
        num_publications[journal_name][str_year_month]['total'] = '25'

        journal_name = 'Infection Genetics and Evolution'
        str_year_month = '2020/1'
        num_publications[journal_name][str_year_month]['total'] = '30'
        # str_year_month = '2020/2'
        # num_publications[journal_name][str_year_month]['total'] = '43'
        str_year_month = '2020/3'
        # num_publications[journal_name][str_year_month]['total'] = '30'

        journal_name = 'Brain  Behavior  and Immunity'
        str_year_month='2020/4'
        num_publications[journal_name][str_year_month]['total'] = '10'

        journal_name = 'Diabetes & Metabolic Syndrome: Clinical Research & Reviews'
        str_year_month = '2020/4'
        num_publications[journal_name][str_year_month]['total'] = '22'

        return num_publications



    def extract_covid_from_scopus(self):
        global data
        offset = 0
        str_offset = str(offset)
        results = []
        for i in range(51):
            data='{"qs":"The",  "pub": "The Lancet",  "loadedAfter": "2020-01-01T00:00:00Z"}'
            # data = '{\n  "pub: \\"The Lancet Infectious Diseases\\" AND 2020",\n "display": {\n      "offset": ' + str_offset + ',\n      "show": 100,\n      "sortBy": "date"\n  }\n}'
            # data = '{\n  "qs": "(\\"COVID-19\\" OR \\"ncov\\") AND 2020",\n "display": {\n      "offset": ' + str_offset + ',\n      "show": 100,\n      "sortBy": "date"\n  }\n}'
            # data = '{\n  "qs": "\\"COVID-19\\" OR \\"ncov\\" AND 2020",\n "display": {\n      "'+str_offset+'": 100,\n      "show": 100,\n      "sortBy": "date"\n  }\n}'
            print(data)
            res = pss.send_request(pss.sciencedirect, pss.headers, data)
            if res!=None:
                results.append(res['results'])
            time.sleep(1)
            offset += 100
            str_offset = str(offset)
        df = self.convert_to_df(results)
        name = 'scopus_updated_3'
        utils.write_to_csv(df, os.path.join(path, name + '.csv'))



if __name__ == '__main__':
    pss=ParseScopusSearch()
    utils=Utils(path=path)
    pss.extract_covid_from_scopus()

    # res_journals = pss.iterate_journals()
    # res_journals = pss.handle_special_journals(res_journals)
    # utils.save_obj(res_journals,"journals_data_apr")
    exit(0)

