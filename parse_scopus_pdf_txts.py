import os
from copy import deepcopy
import numpy as np
from tqdm import tqdm
import re
import pandas as pd


path='D:\\shir\\study\\covid_19\\scopus\\pdfs\\tmp'
rec_in_revised_form='Received in revised form'
rec='Received'
acc='Accepted'
avail_online='Available online'
to_appear_in='To appear in:\n\n'

class PdfTxtsProcess():
    def process_files(self,pdfs_folder_path):
        all_files = []
        print(f"Files in folder: {len(os.listdir(pdfs_folder_path))}")
        # to process all files, uncomment the next line and comment the line below
        # list_of_files = list(os.listdir(xml_folder_path))
        list_of_files = list(os.listdir(pdfs_folder_path))

        venues = []
        dois=[]
        dict_ = {'title': [], 'doi': [], 'venue': [], 'date_received': [], 'date_available':[], 'date_accepted':[] }

        for file in tqdm(list_of_files):
            txt_path = os.path.join(pdfs_folder_path, file)
            if (not txt_path.endswith(".txt")):
                continue

            lines =[]
            try:
                with open(txt_path, encoding="utf8") as txt_file:
                    contents=txt_file.read()

                title = file
                print(title)
                dict_['title'].append(title)
                venue = contents[0:contents.find('\n\n')]
                print(venue)
                if not venue.lower().startswith('journal pre-proof'):
                    if venue.lower().strip()=='articles':
                        venue='Lancet'
                    venue = "".join(filter(lambda x:  x.isalpha() or x.isspace(), venue))
                    dict_['venue'].append(venue)
                else:
                    regex = re.compile('To appear in:\n\n.*', re.IGNORECASE)
                    result = regex.search(contents)
                    if result!=None:
                        venue=contents[result.regs[0][0]+len(to_appear_in):result.regs[0][1]]
                        venue = "".join(filter(lambda x: x.isalpha()or x.isspace(), venue))
                        dict_['venue'].append(venue)
                    else:
                        dict_['venue'].append("unknown")
                regex = re.compile(r'https://doi.org.*', re.IGNORECASE)
                result = regex.search(contents)
                if result!=None:
                    doi = contents[result.regs[0][0]:result.regs[0][1]]
                else:
                    doi='unknown'
                print(doi)
                dois.append(doi)
                dict_['doi'].append(doi)

                regex=re.compile(r'Received \d+ .*; R', re.IGNORECASE)
                result = regex.search(contents)
                if result!=None:
                    date_received=contents[result.regs[0][0]+len(rec):result.regs[0][1]-3]
                    print(date_received)
                else:
                    regex = re.compile(r'Received \d+ .*', re.IGNORECASE)
                    result = regex.search(contents)
                    if result != None:
                        date_received = contents[result.regs[0][0] + len(rec):result.regs[0][1]]
                        print(date_received)
                    else:
                        regex = re.compile(r'Received  Date: \d+ .*', re.IGNORECASE)
                        result = regex.search(contents)
                        if result != None:
                            date_received = contents[result.regs[0][0] + len('Received  Date:'):result.regs[0][1]]
                            print(date_received)
                        else:
                            date_received='unknown'
                date_received = date_received.strip()
                dict_['date_received'].append(date_received)

                regex = re.compile(r'Received in .* \d+ .*; A', re.IGNORECASE)
                result = regex.search(contents)
                if result!=None:
                    date_revised=contents[result.regs[0][0]+len(rec_in_revised_form):result.regs[0][1]-3]
                else:
                    date_revised='unknown'
                print(date_revised)

                regex = re.compile(r'Accepted \d+ .*', re.IGNORECASE)
                result = regex.search(contents)
                if result!=None:
                    date_accepted = contents[result.regs[0][0] + len(acc):result.regs[0][1]]
                else:
                    regex = re.compile(r'Accepted Date: \d+ .*', re.IGNORECASE)
                    result = regex.search(contents)
                    if result != None:
                        date_accepted = contents[result.regs[0][0] + len('Accepted Date:'):result.regs[0][1]]
                    else:
                        date_accepted='unknown'
                date_accepted=date_accepted.strip()
                print(date_accepted)
                dict_['date_accepted'].append(date_accepted)

                regex = re.compile(r'Available online \d+ .*', re.IGNORECASE)
                result = regex.search(contents)
                if result != None:
                    date_available = contents[result.regs[0][0] + len(avail_online):result.regs[0][1]]
                else:
                    date_available = 'unknown'
                date_available=date_available.strip()
                dict_['date_available'].append(date_available)
                print(date_available)

                for venue in venues:
                    print(venue)


            except:
                print("failed for file {}".format(txt_file))

        df_covid = pd.DataFrame(dict_, columns=['title', 'doi', 'venue', 'date_received', 'date_available', 'date_accepted'])
        return df_covid

    def write_to_csv(self, df,name):
        df.to_csv(name, index=False)


if __name__ == '__main__':
    pdf=PdfTxtsProcess()
    df=pdf.process_files(path)
    name='scopus_covid2'
    pdf.write_to_csv(df,os.path.join(path,name+'.csv'))
