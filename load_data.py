import os
import pandas as pd
import json
from copy import deepcopy
import numpy as np
from tqdm import tqdm

class FilesProcess():
    def process_files(self,json_folder_path):
        all_files = []
        print(f"Files in folder: {len(os.listdir(json_folder_path))}")
        # to process all files, uncomment the next line and comment the line below
        # list_of_files = list(os.listdir(json_folder_path))
        list_of_files = list(os.listdir(json_folder_path))


        for file in tqdm(list_of_files):
            json_path = os.path.join(json_folder_path, file)
            with open(json_path) as json_file:
                json_data = json.load(json_file)
            all_files.append(json_data)
            # for key in json_data.keys():
            #     print("key {}, len {}".format(key,len(json_data[key])))
            # json_data_df = pd.DataFrame.from_dict( json_data, orient='index').transpose()
            # df = df.append(json_data_df)
        print("dic keys, ", all_files[0].keys())
        return all_files

    def generate_clean_df(self, all_files):
        cleaned_files = []

        for file in tqdm(all_files):
            features = [
                file['paper_id'],
                file['metadata']['title'],
                self.format_authors(file['metadata']['authors']),
                self.format_authors(file['metadata']['authors'],
                               with_affiliation=True),
                self.format_body(file['abstract']),
                self.format_body(file['body_text']),
                self.format_bib(file['bib_entries']),
                file['metadata']['authors'],
                file['bib_entries']
            ]

            cleaned_files.append(features)

        col_names = ['paper_id', 'title', 'authors',
                     'affiliations', 'abstract', 'text',
                     'bibliography', 'raw_authors', 'raw_bibliography']

        clean_df = pd.DataFrame(cleaned_files, columns=col_names)
        clean_df.head()

        return clean_df

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

    def convert_to_df(self, all_json, meta_df):
        dict_ = {'paper_id': [], 'abstract': [], 'body_text': [], 'authors': [], 'title': [], 'doi': [], 'journal': [],
                 'abstract_summary': []}

        dict_title_doi = {'paper_id': [], 'title': [], 'doi': [], 'authors': []}
        for idx, entry in enumerate(all_json):
            if idx % (len(all_json) // 10) == 0:
                print(f'Processing index: {idx} of {len(all_json)}')
            # content = FileReader(entry)

            # get metadata information
            meta_data = meta_df.loc[meta_df['sha'] == entry['paper_id']]
            # no metadata, skip this paper
            if len(meta_data) == 0:
                continue

            dict_['paper_id'].append(entry['paper_id'])

            abstract = []
            for abs in entry['abstract']:
                abstract.append(abs['text'])
            abstract = '\n'.join(abstract)
            dict_['abstract'].append(abstract)

            body = []
            for body_text in entry['body_text']:
                body.append(body_text['text'])
            body = '\n'.join(body)
            dict_['body_text'].append(body)

            # also create a column for the summary of abstract to be used in a plot

            if len(abstract) == 0:
                # no abstract provided
                dict_['abstract_summary'].append("Not provided.")
            elif len(abstract.split(' ')) > 100:
                # abstract provided is too long for plot, take first 300 words append with ...
                info = abstract.split(' ')[:100]
                summary = self.get_breaks(' '.join(info), 40)
                dict_['abstract_summary'].append(summary + "...")
            else:
                # abstract is short enough
                summary = self.get_breaks(abstract, 40)
                dict_['abstract_summary'].append(summary)

            # get metadata information
            meta_data = meta_df.loc[meta_df['sha'] == entry['paper_id']]

            try:
                # if more than one author
                authors = meta_data['authors'].values[0].split(';')
                if len(authors) > 2:
                    # more than 2 authors, may be problem when plotting, so take first 2 append with ...
                    dict_['authors'].append(". ".join(authors[:2]) + "...")
                else:
                    # authors will fit in plot
                    dict_['authors'].append(". ".join(authors))
            except Exception as e:
                # if only one author - or Null valie
                dict_['authors'].append(meta_data['authors'].values[0])


            # add the title information, add breaks when needed
            try:
                title = self.get_breaks(meta_data['title'].values[0], 40)
                dict_['title'].append(title)

            # if title was not provided
            except Exception as e:
                dict_['title'].append(meta_data['title'].values[0])


            # add the journal information
            dict_['journal'].append(meta_data['journal'].values[0])
            dict_['doi'].append(meta_data['doi'])

        df_covid = pd.DataFrame(dict_, columns=['paper_id', 'abstract', 'body_text', 'authors', 'title','doi', 'journal',
                                                'abstract_summary'])
        df_covid.drop_duplicates(['abstract', 'body_text'], inplace=True)
        df_covid_title_doi= pd.DataFrame()
        df_covid_title_doi['paper_id']=df_covid['paper_id']
        df_covid_title_doi['doi']=df_covid['doi']
        df_covid_title_doi['title']=df_covid['title']
        df_covid_title_doi['authors']=df_covid['authors']



        # df_covid_title_doi.drop_duplicates(['paper_id', 'doi'], inplace=True)
        print(df_covid_title_doi.info())
        return df_covid, df_covid_title_doi

    def write_to_csv(self, df,name):
        df.to_csv(name, index=False)



    def format_name(self, author):
        middle_name = " ".join(author['middle'])

        if author['middle']:
            return " ".join([author['first'], middle_name, author['last']])
        else:
            return " ".join([author['first'], author['last']])

    def format_affiliation(self, affiliation):
        text = []
        location = affiliation.get('location')
        if location:
            text.extend(list(affiliation['location'].values()))

        institution = affiliation.get('institution')
        if institution:
            text = [institution] + text
        return ", ".join(text)

    def format_authors(self, authors, with_affiliation=False):
        name_ls = []

        for author in authors:
            name = self.format_name(author)
            if with_affiliation:
                affiliation = self.format_affiliation(author['affiliation'])
                if affiliation:
                    name_ls.append(f"{name} ({affiliation})")
                else:
                    name_ls.append(name)
            else:
                name_ls.append(name)

        return ", ".join(name_ls)

    def format_body(self, body_text):
        texts = [(di['section'], di['text']) for di in body_text]
        texts_di = {di['section']: "" for di in body_text}

        for section, text in texts:
            texts_di[section] += text

        body = ""

        for section, text in texts_di.items():
            body += section
            body += "\n\n"
            body += text
            body += "\n\n"

        return body

    def format_bib(self, bibs):
        if type(bibs) == dict:
            bibs = list(bibs.values())
        bibs = deepcopy(bibs)
        formatted = []

        for bib in bibs:
            bib['authors'] = self.format_authors(
                bib['authors'],
                with_affiliation=False
            )
            formatted_ls = [str(bib[k]) for k in ['title', 'authors', 'venue', 'year']]
            formatted.append(", ".join(formatted_ls))

        return "; ".join(formatted)




