import requests
import os
import requests
import json
import time
import xmltodict
from utils import Utils
import pycountry
from visualization import Visualization


from Bio import Entrez
max_num_authors=25

data1 = ['10.1016/j.visj.2020.100742', '10.1016/j.visj.2020.100760', '10.1016/S2352-3018(20)30079-5',
         '10.1016/j.knosys.2020.105812', '10.1016/j.ajem.2020.03.052', '10.1016/j.jogc.2020.03.013',
         '10.1016/j.radonc.2020.03.030', '10.1016/j.arbres.2020.03.017', '10.1016/j.adro.2020.03.012',
         '10.1016/j.mjafi.2020.03.009', '10.1016/j.xkme.2020.03.004', '10.1016/j.jtho.2020.03.021',
         '10.1016/j.jogc.2020.03.012', '10.1016/j.jimed.2020.03.001', '10.1016/j.surg.2020.03.012',
         '10.1016/j.jfo.2020.03.001', '10.1016/j.tracli.2020.03.004', '10.1016/j.tifs.2020.03.041',
         '10.1016/j.sapharm.2020.03.020', '10.1016/j.arbres.2020.03.018', '10.1016/j.jaci.2020.03.018',
         '10.1016/j.chmed.2020.03.004', '10.1016/j.clnu.2020.03.022', '10.1016/j.idm.2020.03.004',
         '10.1016/j.jmii.2020.03.022', '10.1016/j.amj.2020.03.003', '10.1016/j.jviscsurg.2020.03.008',
         '10.1016/j.nmni.2020.100672', '10.1016/j.nepr.2020.102780', '10.1016/j.sapharm.2020.03.013',
         '10.1016/j.diii.2020.03.010', '10.1016/j.fsisyn.2020.03.007', '10.1016/j.anpedi.2020.03.002',
         '10.1016/j.jpainsymman.2020.03.018', '10.1016/j.jrid.2020.03.006', '10.1016/j.ijcard.2020.03.074',
         '10.1016/j.tracli.2020.03.005', '10.1016/j.matt.2020.03.017', '10.1016/j.surg.2020.03.011',
         '10.1016/j.jchirv.2020.03.007', '10.1016/j.lpmfor.2020.03.023', '10.1016/j.sheji.2020.02.002',
         '10.1016/j.amp.2020.03.001', '10.1016/j.idm.2020.03.002', '10.1016/S1470-2045(20)30175-3',
         '10.1016/j.jobb.2020.02.001', '10.1016/j.jobb.2020.02.001', '10.1016/j.jobb.2020.02.001',
         '10.1016/j.lpmfor.2020.03.024', '10.1016/j.scitotenv.2020.138201', '10.1016/j.anpedi.2020.03.001',
         '10.1016/j.amcp.2020.02.009', '10.1016/j.amcp.2020.02.009', '10.1016/j.amcp.2020.02.009',
         '10.1016/j.medin.2020.03.001', '10.1016/j.arbres.2020.03.005', '10.1016/j.radonc.2020.03.029',
         '10.1016/j.pharma.2020.03.001', '10.1016/j.wombi.2020.03.008', '10.1016/j.nupar.2020.03.001',
         '10.1016/j.csbj.2020.03.025', '10.1016/j.medin.2020.03.005', '10.1016/j.oftal.2020.02.013',
         '10.1016/j.lpmfor.2020.03.020', '10.1016/j.purol.2020.03.009', '10.1016/j.ijporl.2020.110030',
         '10.1016/j.nupar.2020.03.002', '10.1016/j.lpmfor.2020.03.022', '10.1016/j.kint.2020.03.015',
         '10.1016/j.jphys.2020.03.011', '10.1016/j.healun.2020.03.017', '10.1016/j.healun.2020.03.018',
         '10.1016/j.cjca.2020.03.026', '10.1016/j.cjca.2020.03.028', '10.1053/j.jvca.2020.03.050',
         '10.1016/j.bbih.2020.100064', '10.1016/j.gie.2020.03.3848', '10.1016/j.imr.2020.100407',
         '10.1016/j.eclinm.2020.100325', '10.1016/j.jmii.2020.03.024', '10.1016/j.onehlt.2020.100129',
         '10.1016/j.jmii.2020.03.026', '10.1016/j.resuscitation.2020.03.010', '10.1016/j.jmii.2020.03.021',
         '10.1016/j.bmc.2020.115466', '10.1016/j.ijcard.2020.03.063', '10.1016/j.semerg.2020.03.001',
         '10.1016/j.onehlt.2020.100128', '10.1016/S0262-4079(20)30646-1', '10.1016/S0262-4079(20)30610-2',
         '10.1016/S0262-4079(20)30610-2', '10.1016/S0262-4079(20)30615-1', '10.1016/j.eclinm.2020.100320',
         '10.1016/j.adro.2020.03.007', '10.1053/j.gastro.2020.03.046', '10.1016/j.bulcan.2020.03.001',
         '10.1016/j.bjane.2020.03.002', '10.1016/j.medmal.2020.03.005', '10.1053/j.jvca.2020.03.043',
         '10.1016/S1473-3099(20)30252-8', '10.1016/j.jaccas.2020.03.007', '10.1016/j.msard.2020.102073',
         '10.1016/j.onehlt.2020.100130', '10.1016/j.cdtm.2020.02.003', '10.1016/j.jrid.2020.03.004',
         '10.1016/j.bulcan.2020.03.003', '10.1016/j.medmal.2020.03.004', '10.1016/j.jchf.2020.03.005',
         '10.1016/j.adro.2020.03.006', '10.1016/j.radonc.2020.03.016', '10.1016/j.adro.2020.03.004',
         '10.1016/j.puhip.2020.100004', '10.1016/S1473-3099(20)30182-1', '10.1016/S0140-6736(20)30718-2',
         '10.1016/S0140-6736(20)30719-4', '10.1016/S0140-6736(20)30686-3', '10.1016/S0140-6736(20)30727-3',
         '10.1016/S0140-6736(20)30721-2', '10.1016/j.pedn.2020.03.014', '10.1016/S2213-2600(20)30160-0',
         '10.1016/j.sapharm.2020.03.012', '10.1016/j.patter.2020.100018', '10.1016/j.scitotenv.2020.138226',
         '10.1016/j.xkme.2020.03.003', '10.1016/j.medidd.2020.100028', '10.1016/S2214-109X(20)30115-7',
         '10.1016/j.leukres.2020.106353', '10.1016/j.jasc.2020.03.001', '10.1016/j.jcot.2020.03.014',
         '10.1016/j.aimed.2020.02.003', '10.1016/j.ajogmf.2020.100110', '10.1016/j.ajogmf.2020.100107']
data2 = ['10.1016/j.healun.2020.03.008', '10.1016/j.chb.2020.106357', '10.1016/S1473-3099(20)30251-6',
         '10.1016/j.hrcr.2020.03.012', '10.1016/j.japb.2020.03.012', '10.1016/j.jacbts.2020.03.009',
         '10.1016/S0140-6736(20)30736-4', '10.1016/j.jacr.2020.03.006', '10.1016/j.amcp.2020.03.003',
         '10.1016/j.jmii.2020.03.015', '10.1016/j.gofs.2020.03.017', '10.1016/j.jagp.2020.03.007',
         '10.1016/j.phrs.2020.104768', '10.1016/j.medidd.2020.100026', '10.1016/S2542-5196(20)30062-0',
         '10.1016/j.idm.2020.03.001', '10.1016/j.pdisas.2020.100080', '10.1016/j.tre.2020.101922',
         '10.1016/j.ophtha.2020.03.026', '10.1016/j.clindermatol.2020.03.012', '10.1016/j.ajem.2020.03.036',
         '10.1016/j.shaw.2020.03.001', '10.1016/j.cjca.2020.03.027', '10.1016/j.jacr.2020.03.011',
         '10.1016/j.ctro.2020.03.009', '10.1016/j.buildenv.2020.106827', '10.1016/j.ajem.2020.03.033',
         '10.1016/j.gloepi.2020.100023', '10.1016/j.scitotenv.2020.138149', '10.1016/j.crad.2020.03.003',
         '10.1016/j.jaad.2020.03.037', '10.1016/j.chaos.2020.109761', '10.1016/j.jad.2020.03.041',
         '10.1016/j.jcct.2020.03.002', '10.1016/j.molmed.2020.02.008', '10.1016/j.jmii.2020.03.020',
         '10.1016/j.bjan.2020.03.002', '10.1016/j.resuscitation.2020.03.005', '10.1016/S0262-4079(20)30569-8',
         '10.1016/j.rce.2020.03.001', '10.1016/j.tmaid.2020.101632', '10.1016/j.jhin.2020.03.020',
         '10.1016/j.eng.2020.03.006', '10.1016/j.ajem.2020.03.027', '10.1016/j.jmii.2020.03.019',
         '10.1016/j.healun.2020.03.012', '10.1016/j.hrtlng.2020.03.007', '10.1016/j.adro.2020.03.003',
         '10.1016/j.indmarman.2020.02.017', '10.1016/j.geosus.2020.03.005', '10.1016/S2215-0366(20)30102-4',
         '10.1016/S0140-6736(20)30645-0', '10.1016/S0140-6736(20)30669-3', '10.1016/S0140-6736(20)30672-3',
         '10.1016/S0140-6736(20)30670-X', '10.1016/S0140-6736(20)30644-9', '10.1016/j.ijrobp.2020.03.008',
         '10.1016/j.ijrobp.2020.03.008', '10.1016/j.ijrobp.2020.03.008', '10.1016/j.carrev.2020.03.022',
         '10.1016/j.bbrc.2020.03.047', '10.1016/j.jhin.2020.03.018', '10.1016/j.invent.2020.100317',
         '10.1016/j.ajogmf.2020.100106', '10.1016/j.ijnurstu.2020.103578', '10.1016/j.ijid.2020.03.007',
         '10.1016/j.jtcvs.2020.03.012', '10.1053/j.gastro.2020.03.026', '10.1016/j.scib.2020.03.024',
         '10.1016/j.jmii.2020.03.010', '10.1016/j.fop.2020.03.003', '10.1016/j.jand.2020.02.015',
         '10.1016/S1473-3099(20)30223-1', '10.1016/j.glohj.2020.03.001', '10.1016/j.jaad.2020.03.013',
         '10.1016/S1473-3099(20)30224-3', '10.1016/j.eng.2020.03.002', '10.1016/j.jaad.2020.03.012',
         '10.1016/S0140-6736(20)30629-2', '10.1016/j.wjam.2020.03.005', '10.1016/j.glohj.2020.03.002',
         '10.1016/j.healun.2020.03.006', '10.1016/j.redar.2020.03.003', '10.1016/j.glohj.2020.03.003',
         '10.1016/S2589-7500(20)30059-5', '10.1016/j.ceh.2020.03.001', '10.1016/j.actpha.2020.03.002',
         '10.1016/S2213-2600(20)30128-4', '10.1016/j.jds.2020.02.002', '10.1016/j.dsx.2020.03.005',
         '10.1016/j.fopow.2020.03.001', '10.1016/j.aej.2020.02.033', '10.1016/j.healun.2020.03.007',
         '10.1016/S1473-3099(20)30180-8', '10.1016/S0262-4079(20)30527-3', '10.1016/S0262-4079(20)30522-4',
         '10.1016/S0262-4079(20)30525-X', '10.1016/S0262-4079(20)30547-9', '10.1016/S0262-4079(20)30547-9',
         '10.1016/S0262-4079(20)30547-9', '10.1016/S0262-4079(20)30547-9', '10.1016/S0262-4079(20)30526-1',
         '10.1016/S0262-4079(20)30519-4', '10.1016/S0262-4079(20)30519-4', '10.1016/S0262-4079(20)30519-4',
         '10.1016/S0262-4079(20)30519-4', '10.1016/S0262-4079(20)30519-4', '10.1016/S0262-4079(20)30519-4',
         '10.1016/j.jmii.2020.03.008', '10.1016/j.anpedi.2020.02.001', '10.1016/j.xkme.2020.03.001',
         '10.1016/j.explore.2020.02.022', '10.1016/j.arbres.2020.02.009', '10.1016/S0140-6736(20)30600-0',
         '10.1016/j.jhin.2020.03.010', '10.1016/j.jhin.2020.03.010', '10.1016/j.jhin.2020.03.010',
         '10.1016/j.jhin.2020.03.010', '10.1016/S0966-842X(20)30015-9', '10.1016/j.annonc.2020.03.001',
         '10.1016/j.ajo.2020.02.014', '10.1016/j.cub.2020.02.049', '10.1016/j.jfo.2020.02.001',
         '10.1016/j.jcjq.2020.03.001', '10.1016/j.ctcp.2020.101132', '10.1016/S0262-4079(20)30476-0',
         '10.1016/S0262-4079(20)30472-3', '10.1016/S0262-4079(20)30477-2', '10.1016/S0262-4079(20)30471-1',
         '10.1016/S0262-4079(20)30471-1', '10.1016/S0262-4079(20)30471-1', '10.1016/S0262-4079(20)30471-1',
         '10.1016/S0262-4079(20)30471-1', '10.1016/S0262-4079(20)30471-1', '10.1016/S0262-4079(20)30471-1',
         '10.1016/S0262-4079(20)30471-1', '10.1016/S0262-4079(20)30471-1', '10.1016/S0262-4079(20)30471-1',
         '10.1016/S0262-4079(20)30471-1', '10.1016/S0262-4079(20)30471-1', '10.1016/S0262-4079(20)30474-7',
         '10.1016/j.jpha.2020.03.004', '10.1016/j.jhin.2020.03.001', '10.1016/S0140-6736(20)30522-5',
         '10.1016/S2214-109X(20)30083-8', '10.1016/j.jpha.2020.03.001', '10.1016/j.ctcp.2020.101131',
         '10.1016/j.scib.2020.02.025', '10.1016/j.midw.2020.102668', '10.1016/j.jhin.2020.02.013',
         '10.1016/j.jpha.2020.02.010']


class ParseSpecialAuthors():
    headers = {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'x-els-apikey': 'fee04e1e02e2b7ed9d4991b87f29896f',
    }
    id2='212403'
    id='PMC7104051'



    def get_countries(self, id, countries_dict, countries_rel, doi, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors):
        response = requests.get('https://www.ncbi.nlm.nih.gov//entrez/eutils/efetch.fcgi?db=pmc&id='+id)
        xml_data = xmltodict.parse(response.text)
        try:
            article_meta=xml_data['pmc-articleset']['article']['front']['article-meta']
        except:
            print("from entrez missing part of article {}".format(id))
            unknown.append(doi)
            return countries_dict, countries_rel, doi, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors
        if not 'aff' in article_meta.keys():
            print("from entrez no aff for {}".format(id))
            unknown.append(doi)
            return countries_dict, countries_rel, doi, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors
        affiliations_=article_meta['aff']
        authors=xml_data['pmc-articleset']['article']['front']['article-meta']['contrib-group']
        num_authors=0
        if len(authors)>1:
            num_authors += len(authors)
        else:
            if isinstance(authors['contrib'],list):
                num_authors += len(authors['contrib'])
            else:
                num_authors+=1
        countries=dict()
        if isinstance(affiliations_,list):
            for aff in affiliations_:
                affiliation_=aff['#text']
                # ind_of_dot=affiliation_.find('.')
                ind_of_last_comma=affiliation_.rfind(',',0,-1)
                country=(str(affiliation_)[ind_of_last_comma+1:len(affiliation_)]).strip()

                if not country in countries.keys():
                    countries[country]=0
                countries[country]+=1
        else:
            affiliation_ = affiliations_['#text']
            # ind_of_dot=affiliation_.find('.')
            ind_of_last_comma = affiliation_.rfind(',', 0, -1)
            country = (str(affiliation_)[ind_of_last_comma + 1:len(affiliation_)]).strip()
            countries[country]=num_authors
        if len(countries)==1:
            country_ = self.check_get_country(country)
            if country_ == None:
                country_ = country
            if str(country_).startswith('Taiwan'):
                country_ = 'Taiwan'
            if country_ not in countries_dict:
                countries_dict[country_]=0
                countries_rel[country_]=0
            if country_ not in countries_by_num_articles_dict:
                countries_by_num_articles_dict[country_] = 0
            # print("one country {}".format(country))
            countries_dict[country_]+=num_authors
            countries_rel[country_]+=1
            if num_authors<max_num_authors:
                if num_authors<=10:
                    countries_by_num_authors[num_authors].add(country_)
                else:
                    countries_by_num_authors[num_authors].add(country_+str(doi))
            else:
                countries_by_num_authors[max_num_authors-1].add(country_+str(doi))
            if not country_ in countries_for_article:
                countries_by_num_articles_dict[country_]+=1
                countries_for_article.add(country_)
        else:
            # print("multi countries ")
            for co, val in countries.items():
                country_ = self.check_get_country(co)
                if country_ == None:
                    country_ = co
                if str(country_).startswith('Taiwan'):
                    country_ = 'Taiwan'
                if country_ not in countries_dict:
                    countries_dict[country_] = 0
                    countries_rel[country_] = 0
                if country_ not in countries_by_num_articles_dict:
                    countries_by_num_articles_dict[country_]=0
                countries_dict[country_]+=val
                countries_rel[country_]+=val/num_authors
                if num_authors<max_num_authors:
                    if num_authors<=10:
                        countries_by_num_authors[num_authors].add(country_)
                    else:
                        countries_by_num_authors[num_authors].add(country_+str(doi))
                else:
                    countries_by_num_authors[max_num_authors-1].add(country_+str(doi))
                if not country_ in countries_for_article:
                    countries_by_num_articles_dict[country_]+=1
            handled_current_article=True
        print("success from Entrez")
        return countries_dict, countries_rel, doi, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors


    def get_data_by_doi(self, doi, countries_dict, countries_rel, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors):
            response=requests.get('https://api.elsevier.com/content/article/doi/'+doi+'?apiKey=fee04e1e02e2b7ed9d4991b87f29896f&httpAccept=text%2Fxml')
            response_dict = xmltodict.parse(response.text)
            try:
                try:
                    authors_group=response_dict['full-text-retrieval-response']['originalText']['xocs:doc']['xocs:serial-item']['article']['head']['ce:author-group']
                except:
                    print("failed finding authors for doi {}".format(doi))
                    countries_dict, countries_rel, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors=self.get_data_from_entrez(doi, countries_dict, countries_rel, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors)
                    # unknown.append(doi)
                    return countries_dict, countries_rel, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors
                affilations_dict = dict()
                affilations=authors_group['ce:affiliation']
                authors=authors_group['ce:author']
                num_authors=0
                if isinstance(authors, list):
                    num_authors += len(authors)
                else:
                    num_authors += 1

                # check for single country
                if isinstance(affilations,dict):
                    affi_id = affilations["@id"]
                    country = affilations["sa:affiliation"]["sa:country"]
                    affilations_dict[affi_id]=dict()
                    affilations_dict[affi_id][country] = num_authors
                else:
                    for affi in affilations:
                        affi_id=affi["@id"]
                        if "sa:country" in affi["sa:affiliation"].keys():
                            country=affi["sa:affiliation"]["sa:country"]
                        if "sa:state" in affi["sa:affiliation"].keys():
                            country=affi["sa:affiliation"]["sa:state"]
                        affilations_dict[affi_id] = dict()
                        affilations_dict[affi_id][country]=0
                    if isinstance(authors,list):
                        for author in authors:
                            if isinstance(author['ce:cross-ref'],dict):
                                affil_ref = author['ce:cross-ref']['@refid']
                                for country in affilations_dict[affil_ref].keys():
                                    affilations_dict[affil_ref][country] += 1
                            else:
                                affil_ref = author['ce:cross-ref'][0]['@refid']
                                for country in affilations_dict[affil_ref].keys():
                                    affilations_dict[affil_ref][country]+=1
                    else:
                        affil_ref = authors['ce:cross-ref'][0]['@refid']
                        for country in affilations_dict[affil_ref].keys():
                            affilations_dict[affil_ref][country] += 1
                for k,v in affilations_dict.items():
                    for country, count in v.items():
                        country_=self.check_get_country(country)
                        if country_==None:
                            country_=country
                        if str(country_).startswith('Taiwan'):
                            country_='Taiwan'
                        if country_ not in countries_dict:
                            countries_dict[country_] = 0
                            countries_rel[country_] = 0
                        if country_ not in countries_by_num_articles_dict:
                            countries_by_num_articles_dict[country_]=0
                        if not country_ in countries_for_article:
                            countries_by_num_articles_dict[country_]+=1
                            countries_for_article.add(country_)
                        countries_dict[country_] += count
                        countries_rel[country_] += count/num_authors
                        if num_authors<max_num_authors:
                            if num_authors<=10:
                                countries_by_num_authors[num_authors].add(country_)
                            else:
                                countries_by_num_authors[num_authors].add(country_+str(doi))
                        else:
                            countries_by_num_authors[max_num_authors-1].add(country_+str(doi))
                print('success for doi {}'.format(doi))
            except Exception:
                print('excpetion occured for doi {}'.format(doi))
                countries_dict, countries_rel, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors = self.get_data_from_entrez(doi, countries_dict, countries_rel, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors)
                # unknown.append(doi)
                return countries_dict, countries_rel, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors

            return(countries_dict, countries_rel, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors)

    def get_data_from_entrez(self, doi, countries_dict, countries_rel, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors):
        response = requests.get('https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids='+doi+'&format=json')
        response_dict = json.loads(response.text)
        if 'records' in response_dict.keys() and len(response_dict['records']>0):
            record=response_dict['records'][0]
            if 'pmcid' in record.keys():
                # print(doi)
                pmcid=str(record['pmcid']).strip('PMC')
                countries_dict, countries_rel, doi, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors = self.get_countries(pmcid, countries_dict, countries_rel, doi, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors)
                # print(countries_dict)
                # print(countries_rel)
                # print(response_dict['records'])
            else:
                print("failed from entrez")
                unknown.append(doi)
        else:
            print("failed from entrez")
            unknown.append(doi)
        time.sleep(1)
        return(countries_dict, countries_rel, countries_by_num_articles_dict, unknown, countries_for_article, countries_by_num_authors)


    def check_get_country(self,country):
        try:
            if country.lower().find('hong kong') != -1:
                country = 'Hong Kong'
            if country.lower().find('uk') != -1:
                country = 'United Kingdom'
            if country.lower().find('united states') != -1:
                country = 'USA'
            if country.lower().find('china') != -1:
                country = 'China'
            if country.lower().startswith('espa'):
                country = 'Spain'
            if country.lower().find('singapore') != -1:
                country = 'Singapore'
            if country.lower().find('sydney') != -1:
                country = 'Australia'
            if country.lower().find('taiwan') != -1:
                country = 'Taiwan'
            if country.lower().find('russia') != -1:
                country = 'Russia'
            if country.lower().find('maywood il') != -1:
                country = 'USA'
            if country.lower().find('hackensack') != -1:
                country = 'USA'
            if country.lower().find('maryland') != -1:
                country = 'USA'
            if country.lower().find('alabama') != -1:
                country = 'USA'
            if country.lower().find('albert einstein') != -1:
                country = 'USA'
            if country.lower().find('st louis') != -1:
                country = 'USA'
            if country.lower().find('brunei') != -1:
                country = 'Brunei'
            if country.lower().find('czech') != -1:
                country = 'Czech'
            if country.lower().find('alemania') != -1:
                country = 'Germany'
            if country.lower().find('universitario') != -1:
                country = 'Spain'
            if country.lower().find('new york') != -1:
                country = 'USA'
            if country.lower().find('florida') != -1:
                country = 'USA'
            if country.lower().find('debrecen') != -1:
                country = 'Hungary'
            if country.lower().find('south korea') != -1:
                country = 'Korea'

            res=pycountry.countries.search_fuzzy(country)
            if len(res)>0:
                return res[0].name
            else:
                print("not found {}".format(country))
                return None
        except:
            print("not found {}".format(country))
            return None

    def process_special(self, file_name):
        contents = utils.read_from_csv(file_name).split('\n')
        contents = contents[0:-1]
        countries_dict = dict()
        countries_rel = dict()
        countries_per_article_dict = dict()
        unknowns = []
        num_handled_articles = 0
        for idx, line in enumerate(contents):
            if len(line) == 1:
                continue
            title, authors, doi, journal_name, date, url = line.split(',')
            countries_dict, countries_rel, unknown=self.get_data_by_doi(doi, countries_dict, countries_rel, unknowns)

        # for doi in data1:
        #     countries_dict, countries_rel, unknown=self.get_data_by_doi(doi, countries_dict, countries_rel, unknown)
        # print(len(data1))
        # for doi in data2:
        #     countries_dict, countries_rel, unknown=self.get_data_by_doi(doi, countries_dict, countries_rel, unknown)
        # for key in countries_dict2:
        #     if key in countries_dict1:
        #         countries_dict1[key]+=countries_dict2[key]
        #     else:
        #         countries_dict1[key]=countries_dict2[key]
        # for key in countries_rel2:
        #     if key in countries_rel1:
        #         countries_rel1[key]+=countries_rel2[key]
        #     else:
        #         countries_rel1[key]=countries_rel2[key]
        # unknown=unknown1+unknown2
        # print(len(data2))
        print(len(unknowns))
        print(len(countries_dict))
        a = {k: v for k, v in sorted(countries_dict.items(), key=lambda item: item[1], reverse=True)}
        # for country, count in countries_totals_dict.items():
        #     print(country, count)
        print(a)
        print(unknowns)
        # print("handled articles {}, out of {}".format(num_handled_articles, idx))
        return countries_dict, countries_rel, unknowns


def get_countries_sorted_by_num_articles(file_name):
    global sorted_dict
    res = utils.read_from_csv(file_name).split(',')
    print(res)
    countries_per_title_dict = dict()
    sum = 0
    for idx, line in enumerate(res):
        sep = line.find(':')
        if sep != -1:
            splitted = line.split(':')
            try:
                val = int(splitted[1])
                country = str(splitted[0]).strip('"').strip().strip("'").strip()
                countries_per_title_dict[country] = val
                sum = sum + val
            except:
                print("failed split for {}".format(line))
                continue
        else:
            print(line)
    print(sum)
    sorted_dict = {k: v for k, v in sorted(countries_per_title_dict.items(), key=lambda item: item[1], reverse=True)}
    print(sorted_dict)
    print(idx)
    if 'unknown' in sorted_dict.keys():
        sorted_dict.pop('unknown')
    return sorted_dict


if __name__ == '__main__':
    te=ParseSpecialAuthors()
    utils=Utils(path='D:\\shir\\study\\covid_19\\scopus')
    # file_name='countries per article_Mar.csv'
    # sorted_dict=get_countries_sorted_by_num_articles(file_name)
    # vis=Visualization()
    # vis.plt_country_by_month(sorted_dict, "March")
    # file_name = 'countries per article_Feb.csv'
    # sorted_dict = get_countries_sorted_by_num_articles(file_name)
    # vis = Visualization()
    # vis.plt_country_by_month(sorted_dict, "Feb")
    file_name = 'countries per article_Jan.csv'
    sorted_dict = get_countries_sorted_by_num_articles(file_name)
    vis = Visualization()
    vis.plt_country_by_month(sorted_dict, "Jan")
    file_name = 'countries per article_Feb.csv'
    sorted_dict = get_countries_sorted_by_num_articles(file_name)
    vis = Visualization()
    vis.plt_country_by_month(sorted_dict, "Feb")
    file_name = 'countries per article_Mar.csv'
    sorted_dict = get_countries_sorted_by_num_articles(file_name)
    vis = Visualization()
    vis.plt_country_by_month(sorted_dict, "Mar")
    file_name = 'countries per article_Apr.csv'
    sorted_dict = get_countries_sorted_by_num_articles(file_name)
    vis = Visualization()
    vis.plt_country_by_month(sorted_dict, "Apr")
    # countries_dict1, countries_rel1, unknown=te.process_special('Mar_articles.csv')
    # utils.save_obj(countries_dict1,"pmc_countries_dict_mar")
    # utils.save_obj(countries_rel1,"pmc_countries_rel_mar")
    # utils.save_obj(unknown,"unknown_mar")
    # a=utils.load_obj("unknown_mar")
    # print(len(a))

