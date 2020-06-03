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

max_num_authors=25
class AuthorsAndCountries():
    def __init__(self):
        ## Load configuration
        con_file = open("config.json")
        config = json.load(con_file)
        con_file.close()

        ## Initialize client
        self.client = ElsClient(config['apikey'])
        self.client.inst_token = config['insttoken']
        self.psa = ParseSpecialAuthors()
        self.countries_by_num_authors=[set() for _ in range(max_num_authors)]

    def get_authors_data(self,file_name):

        contents=utils.read_from_csv(file_name).split('\n')
        contents = contents[0:-1]
        authors_dict=dict()
        countries_totals_dict=dict()
        countries_rel_dict=dict()
        countries_per_article_dict=dict()
        handled_current_article=False
        unknowns=[]
        num_handled_articles=0
        for idx, line in enumerate(contents):
            if len(line)==1:
                continue
            title, authors,doi, journal_name,date, url=line.split(',')
            res=self.get_authors_data_by_doi(doi)
            print(res)
            num_authors = len(authors.split(';'))
            if num_authors>1:
                num_authors-=1

            countries_for_article=set()
            if isinstance(res, list):
                num_handled_articles+=1
                for author in res:
                    if len(author['AffiliationInfo'])>0:
                        affiliation_ = author['AffiliationInfo'][0]['Affiliation']
                        ind_of_dot=str(affiliation_).find('.')
                        ind_of_last_comma=str(affiliation_).rfind(',',0,ind_of_dot)
                        country=(str(affiliation_)[ind_of_last_comma+1:ind_of_dot]).strip()
                        if affiliation_.lower().find('hong kong') != -1:
                            country = 'Hong Kong'
                        if country.lower().find('uk') != -1:
                            country = 'United Kingdom'
                        if country.lower().find('united states')!=-1:
                            country='USA'
                        if country.lower().find('china')!=-1:
                            country = 'China'
                        if country.lower().startswith('Espa'):
                            country = 'Spain'
                        if country.lower().find('singapore')!=-1:
                            country = 'Singapore'
                        if country.lower().find('sydney') != -1:
                                country = 'Australia'
                        if country.lower().find('taiwan')!=-1:
                            country = 'Taiwan'
                        country_=self.check_get_country(country)
                        if country_==None:
                            ind_of_next_dot = str(affiliation_).find('.',ind_of_dot+1)
                            ind_of_last_comma = str(affiliation_).rfind(',', ind_of_dot, ind_of_next_dot)
                            country = (str(affiliation_)[ind_of_last_comma + 1:ind_of_next_dot]).strip()
                            if country.lower().startswith('united states'):
                                country = 'USA'
                            if country.lower().startswith('pr china'):
                                country = 'China'
                            if country.lower().startswith('Espa'):
                                country = 'Spain'
                            if country.lower().find('uk') != -1:
                                country = 'United Kingdom'
                            country_ = self.check_get_country(country)
                        if country_==None:
                            if str(affiliation_).lower().find('argentina')!=-1:
                                country='Argentina'
                            if str(affiliation_).lower().find('congo')!=-1:
                                country='Congo'
                            country_ = self.check_get_country(country)
                            if country_==None:
                                country_='unknown'
                        if not country_ in countries_totals_dict.keys():
                            countries_totals_dict[country_]=1
                            countries_rel_dict[country_]=1/len(res)
                            if not country_ in countries_for_article:
                                countries_per_article_dict[country_]=1
                                countries_for_article.add(country_)
                        else:
                            countries_totals_dict[country_]+=1
                            countries_rel_dict[country_]+=1/len(res)
                            if not country_ in countries_for_article:
                                countries_per_article_dict[country_]+=1
                                countries_for_article.add(country_)
                        if len(affiliation_) > 50:
                            affiliation = affiliation_[0:50]
                        if not affiliation in authors_dict.keys():
                            authors_dict[affiliation] = 1
                        else:
                            authors_dict[affiliation] += 1
                        if num_authors<max_num_authors:
                            if num_authors<=10:
                                self.countries_by_num_authors[num_authors].add(country_)
                            else:
                                self.countries_by_num_authors[num_authors].add(country_+str(doi))
                        else:
                            self.countries_by_num_authors[max_num_authors-1].add(country_+str(doi))
                        print('success from pubmed')
                    else:
                        if ((not 'CollectiveName' in author.keys()) or (len(author['CollectiveName']) == 0)):
                            countries_totals_dict, countries_rel_dict, countries_per_article_dict, unknowns, countries_for_article, self.countries_by_num_authors=self.psa.get_data_by_doi(doi, countries_totals_dict, countries_rel_dict, countries_per_article_dict, unknowns, countries_for_article,  self.countries_by_num_authors)

                            # countries_totals_dict, countries_rel_dict = self.handle_specific_doi(doi=doi, countries_totals_dict=countries_totals_dict , countries_rel_dict=countries_rel_dict)
                            # print("unknown author")
                            # unknowns.append("author not known {}".format(doi))

            else:
                countries_totals_dict, countries_rel_dict, countries_per_article_dict, unknowns, countries_for_article, self.countries_by_num_authors=self.psa.get_data_by_doi(doi,countries_totals_dict,countries_rel_dict, countries_per_article_dict, unknowns, countries_for_article, self.countries_by_num_authors)
                # countries_totals_dict, countries_rel_dict = self.handle_specific_doi(doi=doi, countries_totals_dict=countries_totals_dict , countries_rel_dict=countries_rel_dict)
                # print("doi not recognized")
                # unknowns.append("doi not known doi {}".format(doi))

        print(len(countries_totals_dict))
        a={k: v for k, v in sorted(countries_totals_dict.items(), key=lambda item: item[1], reverse=True)}
        # for country, count in countries_totals_dict.items():
        #     print(country, count)
        print(a)
        print(unknowns)
        print("countries per article")
        print(countries_per_article_dict)
        print("countries per author")
        for i,val in enumerate(self.countries_by_num_authors):
            print(i, val)
        # print(*self.countries_by_num_authors, sep=';')
        print("handled articles {}, out of {}".format(num_handled_articles,idx+1))

    def get_authors_data_by_doi(self,doi):
        doi_doc = FullDoc(doi=doi)
        if doi_doc.read(self.client):
            print("doi_doc.title: ", doi_doc.title)
            doi_doc.write()
        else:
            print("Read document failed.")
            return doi
        id=None
        if not 'pubmed-id' in doi_doc._data.keys():
            print("no pubmed-id, trying with title")
            # try with title
            Entrez.email = 'shirAviv@protonmail.com'
            query = doi_doc.title
            handle = Entrez.esearch(db='pubmed',
                                    retmode='xml',
                                    term=query)
            results = Entrez.read(handle)
            if int(results['Count']) > 0:
                id = results['IdList']
        else:
            id = doi_doc._data['pubmed-id']
        if id!=None:
            Entrez.email = 'shirAviv@protonmail.com'
            handle = Entrez.efetch(db='pubmed',
                                   retmode='xml',
                                   id=id)
            results = Entrez.read(handle)
            print(results)
            if 'AuthorList' in results['PubmedArticle'][0]['MedlineCitation']['Article'].keys():
                authors_list = results['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']
                dates = results['PubmedArticle'][0]['PubmedData']['History']
            else:
                print("no authors list {}".format(results['PubmedArticle'][0]['MedlineCitation']['Article']))
                return doi


        else:
            print("no pubmed id")
            return doi

        return authors_list

    def check_get_country(self,country):
        try:
            res=pycountry.countries.search_fuzzy(country)
            if len(res)>0:
                return res[0].name
            else:
                print("not found {}".format(country))
                return None
        except:
            print("not found {}".format(country))
            return None


    def handle_specific_doi_Feb(self, doi, countries_totals_dict, countries_rel_dict):
        country=None
        num_authors=0
        if doi=='10.1016/j.chieco.2020.10143':
            country='USA'
            num_authors=1
        if doi == '10.1016/j.chieco.2020.10143':
            country = 'China'
            num_authors=13
        if doi == '10.1016/j.apsb.2020.02.008':
            country = 'China'
            num_authors=13
        if doi == '10.1016/j.fsiml.2020.100013':
                country = 'Switzerland'
                num_authors = 4
        if doi == '10.1016/j.ajog.2020.02.017':
                country = 'USA'
                num_authors = 4
        if doi == '10.1016/j.tmrv.2020.02.003':
                country = 'China'
                num_authors = 3
        if doi == '10.1016/S2589-7500(20)30026-1':
                country = 'USA'
                num_authors = 3
        if doi == '10.1016/S2213-2600(20)30082-5':
            country = 'United Kingdom'
            num_authors = 1
        if doi == '10.1016/S1473-3099(20)30076-1':
            country = 'United Kingdom'
            num_authors = 1
        if doi == '10.1016/j.scib.2020.02.005':
            country = 'China'
            num_authors = 2
        if doi == '10.1016/j.ceh.2020.02.001':
            country = 'China'
            num_authors = 5
        if doi == '10.1016/j.jshs.2020.01.005':
                country = 'USA'
                num_authors = 1
        if doi == '10.1016/S2213-2600(20)30056-4':
                country = 'United Kingdom'
                num_authors = 1
        country_ = self.check_get_country(country)
        if country_==None:
            country_='unknown'
        if not country_ in countries_totals_dict.keys():
            countries_totals_dict[country_] = num_authors
            if num_authors>0:
                countries_rel_dict[country_] = 1 / num_authors
            else:
                countries_rel_dict[country_] = 0

        else:
            countries_totals_dict[country_] += num_authors
            if num_authors>0:
                countries_rel_dict[country_] += 1 / num_authors

        if doi == '10.1016/j.rce.2020.01.001':
            country = 'USA'
            num_authors = 1
            country_=self.check_get_country(country)
            if country_!=None:
                countries_totals_dict[country_]+=num_authors
                countries_rel_dict[country_] += 1 / num_authors
            else:
                countries_totals_dict['unknown']+=num_authors
                countries_rel_dict['unknown'] += 1 / num_authors

            country='Spain'
            num_authors =1
            country_ = self.check_get_country(country)
            if country_ != None:
                countries_totals_dict[country_] += num_authors
                countries_rel_dict[country_] += 1 / num_authors
            else:
                countries_totals_dict['unknown'] += num_authors
                countries_rel_dict['unknown'] += 1 / num_authors

        if doi == '10.1016/j.hlpt.2020.02.003':
            country = 'United Kingdom'
            num_authors = 1
            country_=self.check_get_country(country)
            if country_!=None:
                countries_totals_dict[country_]+=num_authors
                countries_rel_dict[country_] += 1 / num_authors
            else:
                countries_totals_dict['unknown']+=num_authors
                countries_rel_dict['unknown'] += 1 / num_authors

            country='Netherlands'
            num_authors =1
            country_ = self.check_get_country(country)
            if country_ != None:
                countries_totals_dict[country_] += num_authors
                countries_rel_dict[country_] += 1 / num_authors
            else:
                countries_totals_dict['unknown'] += num_authors
                countries_rel_dict['unknown'] += 1 / num_authors

        return countries_totals_dict, countries_rel_dict

    def handle_specific_doi_Mar(self, doi, countries_totals_dict, countries_rel_dict):
        country=None
        num_authors=0
        if doi=='10.1016/j.visj.2020.100742':
            country='USA'
            num_authors=1
        if doi == '10.1016/j.visj.2020.100760':
            country = 'USA'
            num_authors=2
        if doi == '10.1016/j.knosys.2020.105812':
            country = 'China'
            num_authors=4
        if doi == '10.1016/j.ajem.2020.03.052':
                country = 'USA'
                num_authors = 4

        if doi == '10.1016/j.radonc.2020.03.030':
                country = 'Singapore'
                num_authors = 4
        if doi == '10.1016/j.arbres.2020.03.017':
                country = 'Spain'
                num_authors = 3
        if doi == '10.1016/j.adro.2020.03.012':
                country = 'Switzerland'
                num_authors = 4
        if doi == '10.1016/j.mjafi.2020.03.009':
            country = 'India'
            num_authors = 2
        if doi == '10.1016/j.xkme.2020.03.004':
            country = 'USA'
            num_authors = 2
        if doi == '10.1016/j.jtho.2020.03.021':
            country = 'Italy'
            num_authors = 5

        if doi == '10.1016/j.jimed.2020.03.001':
            country = 'China'
            num_authors = 5
        if doi == '10.1016/j.surg.2020.03.012':
                country = 'Italy'
                num_authors = 6
        if doi == '10.1016/j.jfo.2020.03.001':
                country = 'France'
                num_authors = 2
        if doi == '10.1016/j.jfo.2020.03.001':
                country = 'France'
                num_authors = 2
        if doi == '10.1016/j.jfo.2020.03.001':
            country = 'France'
            num_authors = 2

        if doi == '10.1016/j.tracli.2020.03.004':
                country = 'India'
                num_authors = 1
        if doi == '10.1016/j.tifs.2020.03.041':
                country = 'Italy'
                num_authors = 2
        if doi == '10.1016/j.jfo.2020.03.001':
                country = 'France'
                num_authors = 2
        if doi == '10.1016/j.sapharm.2020.03.020':
                country = 'USA'
                num_authors = 1
        if doi == '10.1016/j.arbres.2020.03.018':
                country = 'Spain'
                num_authors = 3
        if doi == '10.1016/j.jaci.2020.03.018':
                country = 'USA'
                num_authors = 4
        if doi == '10.1016/j.chmed.2020.03.004':
                country = 'China'
                num_authors = 1
        if doi == '10.1016/j.clnu.2020.03.022':
                country = 'France'
                num_authors = 2


        country_ = self.check_get_country(country)
        if country_==None:
            country_='unknown'
        if not country_ in countries_totals_dict.keys():
            countries_totals_dict[country_] = num_authors
            if num_authors>0:
                countries_rel_dict[country_] = 1 / num_authors
            else:
                countries_rel_dict[country_] = 0

        else:
            countries_totals_dict[country_] += num_authors
            if num_authors>0:
                countries_rel_dict[country_] += 1 / num_authors

        if doi == '10.1016/j.clnu.2020.03.022':
            country = 'Germany'
            num_authors = 2
            country_=self.check_get_country(country)
            if country_!=None:
                countries_totals_dict[country_]+=num_authors
                countries_rel_dict[country_] += 1 / num_authors
            else:
                countries_totals_dict['unknown']+=num_authors
                countries_rel_dict['unknown'] += 1 / num_authors

            country='Italy'
            num_authors =1
            country_ = self.check_get_country(country)
            if country_ != None:
                countries_totals_dict[country_] += num_authors
                countries_rel_dict[country_] += 1 / num_authors
            else:
                countries_totals_dict['unknown'] += num_authors
                countries_rel_dict['unknown'] += 1 / num_authors

            country = 'Croatia'
            num_authors = 1
            country_ = self.check_get_country(country)
            if country_ != None:
                countries_totals_dict[country_] += num_authors
                countries_rel_dict[country_] += 1 / num_authors
            else:
                countries_totals_dict['unknown'] += num_authors
                countries_rel_dict['unknown'] += 1 / num_authors

            country = 'Israel'
            num_authors = 1
            country_ = self.check_get_country(country)
            if country_ != None:
                countries_totals_dict[country_] += num_authors
                countries_rel_dict[country_] += 1 / num_authors
            else:
                countries_totals_dict['unknown'] += num_authors
                countries_rel_dict['unknown'] += 1 / num_authors

        if doi == '10.1016/j.hlpt.2020.02.003':
            country = 'United Kingdom'
            num_authors = 1
            country_=self.check_get_country(country)
            if country_!=None:
                countries_totals_dict[country_]+=num_authors
                countries_rel_dict[country_] += 1 / num_authors
            else:
                countries_totals_dict['unknown']+=num_authors
                countries_rel_dict['unknown'] += 1 / num_authors

            country='Netherlands'
            num_authors =1
            country_ = self.check_get_country(country)
            if country_ != None:
                countries_totals_dict[country_] += num_authors
                countries_rel_dict[country_] += 1 / num_authors
            else:
                countries_totals_dict['unknown'] += num_authors
                countries_rel_dict['unknown'] += 1 / num_authors

        return countries_totals_dict, countries_rel_dict


if __name__ == '__main__':
    print(datetime.datetime.now())
    utils=Utils(path='D:\\shir\\study\\covid_19\\scopus')
    pa=AuthorsAndCountries()
    pa.get_authors_data('Apr_articles.csv')
    print(datetime.datetime.now())


