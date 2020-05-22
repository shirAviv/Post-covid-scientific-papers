import time
import numpy as np
import pandas as pd
import os

class Metadata():
    def load_metadata(self, path):
        t1 = time.time()
        df = pd.read_csv(path,dtype={'pubmed_id': str, 'Microsoft Academic Paper ID': str, 'doi': str })
        t2 = time.time()
        print('Elapsed time:', t2 - t1)
        df.drop_duplicates(subset=['title'], inplace=True)
        df.doi = df.doi.fillna('').apply(self.doi_url)
        df.info()
        return df

    def doi_url(self,d):
        if d.startswith('http://'):
            return d
        elif d.startswith('doi.org'):
            return f'http://{d}'
        else:
            return f'http://doi.org/{d}'

    def get_general_stats(self, df):
        print(df.head())
        print(df.describe(include='all'))
        print(df.journal.value_counts())
        print(df.source_x.value_counts())


    def get_publish_year_stats(self, df, venue='all', filter_by=None):
        if venue=='all':
            venue_data=df.publish_time
        else:
            venue_data=df['source_x']==venue
            venue_data=df[venue_data].publish_time
        if filter_by!=None:
            a = venue_data[venue_data.str.count(r'(^'+filter_by+'.*)') == 1]
            print(len(a))
        else:
            a= venue_data
            print(len(a))
        return venue_data
