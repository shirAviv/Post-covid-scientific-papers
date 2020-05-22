
import numpy as np
from pprint import pprint

import os

class Stats():

    def top_level_stats(self, path):
        count = 0
        file_exts = []
        for dirname, _, filenames in os.walk(path):
            for filename in filenames:
                count += 1
                file_ext = filename.split(".")[-1]
                file_exts.append(file_ext)

        file_ext_set = set(file_exts)

        print(f"Files: {count}")
        print(f"Files extensions: {file_ext_set}\n\n=====================\nFiles extension count:\n=====================")
        file_ext_list = list(file_ext_set)
        for fe in file_ext_list:
            fe_count = file_exts.count(fe)
            print(f"{fe}: {fe_count}")

        count = 0
        for root, folders, filenames in os.walk(path):
            print(root, folders)

    def get_venues(self,all_files):
        venues=dict()
        for file in all_files:
        #     pprint(file['metadata']['title'])
        #     pprint('...bib venues....')
            bibs = list(file['bib_entries'].values())
            for bib in bibs:
                venue=str(bib['venue']).replace(' ','_').replace('.','').lower()
                if venue in venues.keys():
                    venues[venue]+=1
                else:
                    venues[venue]=1

        for key, venue in venues.items():
            print(f"{key}: {venue}")
        return venues

    # def get_publication_year(self,all_files):
    #     pub_years=dict()
    #     for file in all_files:

