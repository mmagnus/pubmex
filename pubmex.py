#! /usr/bin/env python

"""
copyright 2010 Marcin Magnus

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

from Bio import Entrez#biopython
import sys
import optparse
import os
import re

VERSION = '0.02'
MAIL='your_mail@gmail.com'

def dot(text):
    """ change ' '(space) into . ,e.g. Marcin Magnus -> Marcin.Magnus"""
    return text.replace(' ','.')

def clean_string(text):
    """
    execute operations to clean up a given string
    """
    text = text.replace(')','.')
    text = text.replace('(','.')
    text = text.replace('..','.') # Marcin..Magnus -> Marcin.Magnus # to do
    text = text.replace('-.','-') # Magnus-.student -> Magnus-student
    text = text.replace('.-','-')
    return text

def prepare_customed_title(text):
    """
    GET:
    - text = customed_title, e.g. 'RNA, structure' or 'RNA structure'
    """
    text = text.replace(',',' ')
    text = re.sub('\s+','.',text)
    return text 

def is_it_pmid(id):
    """
    id - it might be doi or pmid
    """
    if re.search('^\d+$',id):#17123955
        return True 
    else:
        return False
    

def get_title_via_pmid(pmid, reference, customed_title, verbose = 0):
    """
    GET:
    - pmid of a publication, e.g. 17123955
    and so on..
    - v(erbose) 
    DO:
    - use biopython to get summary dict 
    
    summary_dict 
    
    {'DOI': '10.1261/rna.283207', 'Title': 'Conversion of pre-RISC to holo-RISC by Ago2 during assembly of RNAi complexes.', 'Source': 'RNA', 'PmcRefCount': 9, 'Issue': '1', 'SO': '2007 Jan;13(1):22-9', 'ISSN': '1355-8382', 'Volume': '13', 'FullJournalName': 'RNA (New York, N.Y.)', 'RecordStatus': 'PubMed - indexed for MEDLINE', 'ESSN': '1469-9001', 'ELocationID': '', 'Pages': '22-9', 'PubStatus': 'ppublish+epublish', 'AuthorList': ['Kim K', 'Lee YS', 'Carthew RW'], 'EPubDate': '2006 Nov 22', 'PubDate': '2007 Jan', 'NlmUniqueID': '9509184', 'LastAuthor': 'Carthew RW', 'ArticleIds': {'pii': 'rna.283207', 'medline': [], 'pubmed': ['17123955'], 'pmc': 'PMC1705758', 'pmcid': 'pmc-id: PMC1705758;', 'doi': '10.1261/rna.283207'}, u'Item': [], 'History': {'medline': ['2007/01/31 09:00'], 'pubmed': ['2006/11/25 09:00'], 'entrez': '2006/11/25 09:00', 'aheadofprint': '2006/11/22 00:00'}, 'LangList': ['English'], 'HasAbstract': 1, 'References': [], 'PubTypeList': ['Journal Article'], u'Id': '17123955'}

    RETURN:
    - formated title(string) 

    need xclip - http://sourceforge.net/projects/xclip/
    """
    result = Entrez.esummary(db = "pubmed", id = pmid, email = MAIL)
    summary_dict = Entrez.read(result)[0]    

    if verbose:
        print summary_dict

    if not reference:
        title_of_pub = dot(summary_dict['Title'].strip())

        first_author = summary_dict['AuthorList'][0].split(' ')[0].strip()
        title =  first_author

        last_author = summary_dict['LastAuthor'].split(' ')[0].strip()
        if first_author != last_author:
            title += '.'+ last_author

        if not customed_title:
            title += '-'+title_of_pub
        else:
            title += '-'+prepare_customed_title(customed_title)
        title += '-'+pmid
        title += '.'+dot(summary_dict['Source']).upper()
        title += '.'+str(summary_dict['PubDate'].split()[0])
        title += '.pdf'

        title = clean_string(title)
    else:##reference
        title = ", ".join(summary_dict['AuthorList']) + " " + summary_dict['Title'].strip() + " "\
            + summary_dict['Source'].strip() +" "\
            + summary_dict['SO']

    cmd = "echo '"+title + "' | xclip -selection clipboard"
    if verbose:
        print cmd
    os.system(cmd)

    return title

def get_title_via_doi(doi, reference, customed_title, verbose = 0):
    """
    GET:
    - doi, e.g.
    and so on..
    - v(erbose) 
    DO:
    - search in pubmed a paper based on given doi,
    - get the paper, fetch pmid
    - get_title_via_pmid
    RETURN:
    - title(string)
    """
    if verbose:
        print doi

    result = Entrez.esearch(db = "pubmed", term = doi, email = MAIL) ## *need to be improved*
    out = Entrez.read(result)
    idlist = out["IdList"]
    return get_title_via_pmid(idlist[0],reference, customed_title)

def get_options(verbose=False):
    """
    get options
    """
    usage = "%prog [options] id"
    description = """
examples: pubmex.py -p 17123955, pumex.py -d 10.1038/embor.2008.212
"""
    version = VERSION

    parser = optparse.OptionParser(description = description,
                              version = version,
                              usage = usage)

    parser.add_option("--pubmed_id", "-p",
                      help = "pass PMID of the paper")

    parser.add_option("--customed_title", "-t", type="string", #dest='customed_title', 
                      help = """pass your title for a pdf
in a format 'RNA, structure' """ )

    parser.add_option("--doi", "-d",
                      help = "pass DOI of the paper")

    parser.add_option("--reference", "-r", action="store_true",
                      help = "reference format", default=False)
    if verbose:
        print parser.parse_args()#(<Values at 0xb7b4b4cc: {'pmid': '1212'}>, [])

    (options, arguments) = parser.parse_args()

    if not options.pubmed_id and not options.doi:
        parser.print_help()
        parser.error('You have to pass PMID or DOI')
        sys.exit(1)

    return options, arguments

if '__main__' == __name__:
    OPTIONS, ARGUMENTS = get_options()
    if OPTIONS.pubmed_id:
        #if ARGUMENTS[0] == 'demo':
        #    ARGUMENTS[0] = '17123955'
        #    print 'Demo PMID: ',ARGUMENTS[0]
        print get_title_via_pmid(OPTIONS.pubmed_id, OPTIONS.reference, OPTIONS.customed_title)
    if OPTIONS.doi:
        print get_title_via_doi(OPTIONS.doi, OPTIONS.reference, OPTIONS.customed_title)
