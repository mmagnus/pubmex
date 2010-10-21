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

DEP:
 - pdftotexthttp://en.wikipedia.org/wiki/Pdftotext
 - xclip http://sourceforge.net/projects/xclip/
"""

from Bio import Entrez#biopython
import sys
import optparse
import os
import re
import tempfile
import subprocess
import shutil

VERSION = '0.02'
MAIL='your_mail@gmail.com'
DEBUG = False
JDICT={'NUCLEIC.ACIDS.RES': 'NAR'}

def dot(text):
    """ change ' '(space) into . ,e.g. Marcin Magnus -> Marcin.Magnus"""
    return text.replace(' ','.')

def clean_string(text):
    """
    execute operations to clean up a given string
    """
    text = text.replace(')','.')
    text = text.replace('[','.')
    text = text.replace(']','.')
    text = text.replace(':.','.')# for windows, there can not be :
    text = text.replace('(','.')
    text = text.replace('-.','-') # Magnus-.student -> Magnus-student
    text = text.replace('.-','-')
    text = text.replace(',.','.')
    text = text.replace('/','.')
    text = text.replace('\'','')
    text = text.replace('--','-')
    text = text.replace('..','.') # Marcin..Magnus -> Marcin.Magnus # to do
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
    try:
        summary_dict = Entrez.read(result)[0]    
    except RuntimeError: #### there is NO pmd like this
        return False

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
        ###
        for k in JDICT.keys():
            title = title.replace(k, JDICT[k])
            
    else:##reference
        title = ", ".join(summary_dict['AuthorList']) + " " + summary_dict['Title'].strip() + " "\
            + summary_dict['Source'].strip() +" "\
            + summary_dict['SO']


    return title

def text2clip(title, verbose = 0):
    cmd = "echo '"+title.strip() + "'| xclip -selection clipboard"
    if verbose:
        print cmd
    #print 'The filename went to your clipboard :-) (hopefully)'
    os.system(cmd)

def doi2pmid(doi):
    """
    Get a DOI based on a give DOI
    """
    result = Entrez.esearch(db = "pubmed", term = doi, email = MAIL) ## *need to be improved*
    out = Entrez.read(result)
    idlist = out["IdList"]
    if DEBUG: print 'IdList', idlist
    if len(idlist) == 1:
        return idlist[0] ## if IdList ['20959296', '20959295', '20959294', '20959293', '20959292', '20959291', '20959290', '20959289', '20959288', '20959287', '20952411', '20952410', '20952409', '20952408', '20952407', '20952406', '20952405', '20952404', '20952403', '20952402'] # problem !!!
    else: 
        return False

def get_title_auto_from_pdf(filename, reference, customed_title, verbose = 0):
    doi = get_doi_from_pdf(filename)
    if doi:
        return get_title_via_doi(doi, reference, customed_title, verbose = 0)
    else:
        print 'DOI has *not* been found automatically!'
        return False

def get_doi_from_pdf(filename, verbose = 1):
    f = tempfile.NamedTemporaryFile()
    #print f.name
    ## degum
    if not DEBUG:
        txtfn = f.name
    else:
        txtfn = 'temp'
        print 'DEBUG - temp generated'
    f.close()
    #args = ['pdftotext', filename, txtfn]#f.name]
    args = ['pdftotext', filename, txtfn]
    p = subprocess.call(args)
    doi = False
    if p == 0:# it means it's OK
        txt = open(txtfn).read()
        ###
        rex = re.compile('doi:(?P<doi>.*\.\w+\.\w+)').search(txt)
        rex2 = re.compile('DOI\s+(?P<doi>.*\.\w+\.\w+)').search(txt)
        rex3 = re.compile('DOI:\s+(?P<doi>.*\.\w+\.\w+)').search(txt)
        rex4 = re.compile('doi:(?P<doi>\d+.\d+/[\w/]+)').search(txt)
        rex5 = re.compile('DOI:\s+(?P<doi>.+)').search(txt)
        if rex: doi = rex.group('doi')
        if rex2: doi = rex2.group('doi')
        if rex3: doi = rex3.group('doi')
        if rex4: doi = rex4.group('doi')
        if rex5: doi = rex5.group('doi')
        if verbose: print 'DOI:', doi

        ###
    return doi


def get_title_via_doi(doi, reference, customed_title, verbose = 1):
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
    doi=doi.replace('DOI:','')
    doi=doi.replace('doi:','')
    if DEBUG: print doi

    pmid = doi2pmid(doi)
    if DEBUG: print 'PMID:', pmid
    if pmid:
        return get_title_via_pmid(pmid,reference, customed_title)
    else: ### return problem
        return False
        

def get_options(verbose=False):
    """
    get options
    """
    usage = "%prog [options] id"
    description = """
examples: pubmex.py -p 17123955; pumex.py -p 10.1038/embor.2008.212; pubmex.py -a -f file.pdf -r
"""
    version = VERSION

    parser = optparse.OptionParser(description = description,
                              version = version,
                              usage = usage)

    parser.add_option("--PMID_or_DOI", "-p",
                      help = "pass PMID/DOI of the paper")

    parser.add_option("--filename", "-f",
                      help = "filename of a pdf")

    parser.add_option("--automatic", "-a", action="store_true",
                      help = "try to get DOI automatically from a pdf")

    parser.add_option("--keywords", "-k", type="string", #dest='customed_title', 
                      help = """pass your keywords which makes a filename in a format 'RNA, structure' """ )

    #parser.add_option("--reference", "-r", action="store_true",
    #                  help = "reference format", default=False)
    
    parser.add_option("--rename", "-r", action="store_true",
                      help = "rename files (only in a automatic mode)", default=False)


    if verbose:
        print parser.parse_args()#(<Values at 0xb7b4b4cc: {'pmid': '1212'}>, [])

    (options, arguments) = parser.parse_args()
    if DEBUG: print options

    if not options.PMID_or_DOI and not options.automatic:
        parser.print_help()
        parser.error('You have to pass PMID or DOI. You migth set -a (automatic)')
        sys.exit(1)

    return options, arguments

def rename(src, dst, rename_flag):
            if hashverb: print '###'
            if rename_flag:
                print src,'-*DO*->', dst
                shutil.move(src, dst)
            else:
                print src,'-*NOT*->', dst
                pass

if '__main__' == __name__:
    ###
    hashverb = False
    ###
    OPTIONS, ARGUMENTS = get_options()
    keywords = ''
    title = ''
    print '-' * 40
    ###
    if hashverb: print '###'
    if OPTIONS.PMID_or_DOI:
        if is_it_pmid(OPTIONS.PMID_or_DOI):
            title = get_title_via_pmid(OPTIONS.PMID_or_DOI, False, OPTIONS.keywords)
        else:

            title = get_title_via_doi(OPTIONS.PMID_or_DOI, False, OPTIONS.keywords)# OPTIONS.reference, 
    ###
    if OPTIONS.automatic:
        filename = OPTIONS.filename
        print 'FILENAME:', 
        title = get_title_auto_from_pdf(filename, False, OPTIONS.keywords)
        if title:
            basename =  os.path.basename(filename)
            src = filename
            dst = filename.replace(basename, title)
            rename(src, dst, OPTIONS.rename)
        else:
            ## title = False
            pass
    ###
    if title is not False:
        print title
        if hashverb: print '###'
        try:
            text2clip(title)
        except:
            pass
        if OPTIONS.filename and not OPTIONS.automatic:
            basename =  os.path.basename(OPTIONS.filename)
            src = OPTIONS.filename
            dst = OPTIONS.filename.replace(basename, title)
            rename(src, dst, OPTIONS.rename)
    else:
        print 'Problem! Check your PMID/DOI'
        if hashverb: print '###'
