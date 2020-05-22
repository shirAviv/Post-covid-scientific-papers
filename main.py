import os
from load_data import FilesProcess
from get_stats import Stats
from metadata_evaluation import Metadata
from visualization import Visualization
from pprint import pprint

path='D:\\shir\\study\\covid_19\\CORD-19-research-challenge\\2020-03-13'

if __name__ == '__main__':
    metadata=Metadata()
    md=metadata.load_metadata(path+'\\all_sources_metadata_2020-03-13.csv')

    # metadata.get_general_stats(md)
    # venue_data= metadata.get_publish_year_stats(md,venue='PMC',filter_by='2020')
    # visualize=Visualization()
    # visualize.show_year(md)
    # visualize.show_sources(md)
    # exit(0)
    stats=Stats()
    stats.top_level_stats(path)
    data=FilesProcess()
    name='biorxiv_medrxiv'
    all_files=data.process_files(os.path.join(path,name,name))
    df, df_title_doi=data.convert_to_df(all_files,md)
    # df=data.generate_clean_df(all_files)
    data.write_to_csv(df,os.path.join(path,name+'.csv'))
    data.write_to_csv(df_title_doi,os.path.join(path,name+'_title_doi.csv'))

    stats.get_venues(all_files)
    # pprint('...abstracts....')
    # for file in all_files:
    #     pprint(file['abstract'])

    # pprint('...title....')
    # for file in all_files:
    #     pprint(file['metadata']['title'])
    #     pprint('...bib venues....')
    #     bibs = list(file['bib_entries'].values())
    #     for bib in bibs:
    #         pprint(bib['venue'])
