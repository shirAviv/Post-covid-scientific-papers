import os
import pandas as pd
import xml
import xmltodict
from copy import deepcopy
import numpy as np
from tqdm import tqdm

path='D:\\shir\\study\\covid_19\\scopus\\xmls'

class XmlProcess():
    def process_files(self,xml_folder_path):
        all_files = []
        print(f"Files in folder: {len(os.listdir(xml_folder_path))}")
        # to process all files, uncomment the next line and comment the line below
        # list_of_files = list(os.listdir(xml_folder_path))
        list_of_files = list(os.listdir(xml_folder_path))


        for file in tqdm(list_of_files):
            xml_path = os.path.join(xml_folder_path, file)

            try:
                with open(xml_path) as xml_file:
                    xml_data = xmltodict.parse(xml_file.read())
                print("success")

                # y = BeautifulSoup(xml_data)

                # all_files.append(xml_data)
                orig_load_date=xml_data["full-text-retrieval-response"]["originalText"]["xocs:doc"]["xocs:meta"]["xocs:orig-load-date"]['#text']
                if orig_load_date.startswith('2020'):
                    available_online_date=xml_data["full-text-retrieval-response"]["originalText"]["xocs:doc"]["xocs:meta"]["xocs:available-online-date"]['#text']
                    title=xml_data["full-text-retrieval-response"]["coredata"]["dc:title"]
                    publication_name=xml_data["full-text-retrieval-response"]["originalText"]["xocs:doc"]["xocs:meta"]["xocs:srctitle"]
                    doi=xml_data["full-text-retrieval-response"]["originalText"]["xocs:doc"]["xocs:meta"]["xocs:doi"]
                    file_name=file
                    print("found, name - {}, orig date{}, available on line date {}, title {}, venue {}".format(file_name,orig_load_date,available_online_date,title,publication_name))
            except:
                print("skipping file {}".format(file))


if __name__ == '__main__':
    xml_p=XmlProcess()
    xml_p.process_files(path)

