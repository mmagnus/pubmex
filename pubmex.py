#! /usr/bin/env python
#-*- coding: utf-8 -*-
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
 - pdftotext http://en.wikipedia.org/wiki/Pdftotext
 - xclip http://sourceforge.net/projects/xclip/
 - biopython http://biopython.org/wiki/Biopython
"""

__version__ = '0.3a'

from Bio import Entrez
import sys
import optparse
import os
import re
import tempfile
import subprocess
import shutil
import urllib

MAIL = 'your_mail@gmail.com'
DEBUG = True
JDICT = {'NUCLEIC.ACIDS.RES': 'NAR'}
ADD_PMID = False
WORDS_TO_REMOVE = 'a, as, at, for, from, he, her, his, if, in, it, its, of, on, she, so, the, their, them, they, to, which' + ',with, and, by, during'
DONE_MOVE_TO_FOLDER = True
DONE_FOLDER_NAME = 'done'
###############################
if not DEBUG:
    f = tempfile.NamedTemporaryFile()
    TEMPFILE_NAME = f.name
else:
    TEMPFILE_NAME = 'temp'
    print 'DEBUG - temp generated'
###############################


def hr():
    """Draw -------------."""
    print '-' * 80


def dot(text):
    """Change ' '(space) into . ,e.g. Marcin Magnus -> Marcin.Magnus."""
    return text.replace(' ', '.')


def clean_string(text):
    """Execute operations to clean up a given string."""
    for word in WORDS_TO_REMOVE.split(','):
        text = text.replace('.' + word.strip() + '.', '.')
        text = text.replace('.' + word.strip().upper() + '.', '.')
    text = text.replace(')', '.')
    text = text.replace('[', '.')
    text = text.replace(']', '.')
    text = text.replace(':.', '.')  # for windows, there can not be :
    text = text.replace('(', '.')
    text = text.replace('-.', '-')  # Magnus-.student -> Magnus-student
    text = text.replace('.-', '-')
    text = text.replace(',.', '.')
    text = text.replace('/', '.')
    text = text.replace('\'', '')
    text = text.replace('--', '-')
    text = text.replace('..', '.')  # Marcin..Magnus -> Marcin.Magnus # to do
    return text


def prepare_customed_title(text):
    """Only take text e.g 'RNA, structure' and return 'RNA.structure'

    input:
     * text: customed_title
    """
    text = text.replace(',', ' ')
    text = re.sub('\s+', '.', text)
    return text


def is_it_pmid(id):
    """Check if id is PMID."""
    if re.search('^\d+$', id):  # 17123955
        return True
    else:
        return False


#problem?!
def get_pmid_via_search_in_pubmex_line_by_line(
    text='RNA tertiary structure prediction with ModeRNA.'):
    """
    http://baoilleach.blogspot.com/2008/02/searching-pubmed-with-python.html
    """
    ### branie tych lini chyba powinno byc w innej lini !!!!
    #funkcja ktora zczytuje kolejen linie a inna aby
    # odczytywac czy dla naej lini jest unikalny wpis w pubmed
    #1. dla lini pubmed
    #2. dla unikalnego pmid sprawdz czy nazwiska autorow
    # sa w 10 pierwszych liniach textu! !!!
    #txtfn = TEMPFILE_NAME
    #args = ['pdftotext', filename, txtfn]
    #p = subprocess.call(args)
    #doi = False
    #if p == 0:# it means it's OK
    #    txt = open(txtfn).read()
    #    if verbose: print txtfn, "is going to be opend"
    result = Entrez.esearch(db="pubmed", term=text, email=MAIL)
    d = Entrez.read(result)
    print 'd[Count]: ', d['Count']
    pmid = d['IdList']
    print 'pmid: ', pmid
    result = Entrez.esummary(db="pubmed", id=pmid, email=MAIL)
    summary_dict = Entrez.read(result)[0]
    print 'summary_dict: ', summary_dict
    print 'summary_dict[AuthorList]: ', summary_dict['AuthorList']


def get_title_via_pmid(pmid, reference, customed_title, verbose=0):
    """Use biopython to get summary dict.

    input:
     * pmid of a publication, e.g. 17123955
       and so on..
     * v(erbose)

    output:
     * formated title(string)
    """
    result = Entrez.esummary(db="pubmed", id=pmid, email=MAIL)
    try:
        summary_dict = Entrez.read(result)[0]
    except RuntimeError:  # there is NO pmd like this
        return False

    if verbose:
        print 'summary_dict', summary_dict

    if not reference:
        title_of_pub = dot(summary_dict['Title'].strip())

        first_author = summary_dict['AuthorList'][0].split(' ')[0].strip()
        title = first_author

        last_author = summary_dict['LastAuthor'].split(' ')[0].strip()
        if first_author != last_author:
            title += '.' + last_author

        if not customed_title:
            title += '-' + title_of_pub
        else:
            title += '-' + prepare_customed_title(customed_title)

        if ADD_PMID:
            title += '-' + pmid

        title += '-' + dot(summary_dict['Source'].replace(' ', ''))
        title += '.' + str(summary_dict['PubDate'].split()[0])
        title += '.pdf'

        title = clean_string(title)

        for k in JDICT.keys():
            title = title.replace(k, JDICT[k])
    # TODO refernece
    else:
        title = ", ".join(summary_dict['AuthorList']) + " " + summary_dict['Title'].strip() + " " + summary_dict['Source'].strip() + " " + summary_dict['SO']
    return title


def text2clip(title, verbose=0):
    """Output is sent to the clipboard."""
    cmd = "echo '" + title.strip() + "'| xclip -selection clipboard"
    os.system(cmd)


def doi2pmid(doi):
    """Get a PMID based on a given DOI."""
    # *needto be improved*
    result = Entrez.esearch(db="pubmed", term=doi, email=MAIL)
    out = Entrez.read(result)
    idlist = out["IdList"]
    if DEBUG:
        print 'IdList', idlist
    if len(idlist) == 1:
        return idlist[0]
        # if IdList ['20959296', '20959295', '20959294'959289', '20959288',
        #'20959287', '20952411', '20952403', '20952402'] # problem !!!
    else:
        return False


def get_title_auto_from_text(text, reference, customed_title, verbose=0):
    doi = get_doi_from_text(text)
    if doi:
        return get_title_via_doi(doi, reference, customed_title, verbose=0)
    else:
        print 'DOI has *not* been found automatically!'
        return False

def pdf2text(filename):
    """Convert a pdf file to a flat file.

    input:
    * filename

    caution:
    * pdftotext: is required
    """
    txtfn = TEMPFILE_NAME
    args = ['pdftotext', filename, txtfn]
    p = subprocess.call(args)
    doi = False
    if p == 0:
        txt = open(txtfn).read()
        if DEBUG:
            print txtfn, "is going to be opened"
    return txt
    
def get_doi_from_text(text, verbose=False):
    """Use several regular expression are used to get DOI in a text.

    input:
    * text
    * verbose

    output:
    * doi
    """
    ### TODO
    #rex = re.compile('doi:(?P<doi>.*\.\w+\.\w+)').search(txt)
    #rex2 = re.compile('DOI\s+(?P<doi>.*\.\w+\.\w+)').search(txt)
    #rex3 = re.compile('DOI:\s+(?P<doi>.*\.\w+\.\w+)').search(txt)
    #rex4 = re.compile('doi:(?P<doi>\d+.\d+/[\w/]+)').search(txt)
    #rex5 = re.compile('DOI:\s+(?P<doi>.+)').search(txt)
    ### cleaning up the text
    text = text.upper()
    text = text.replace('–', '-')
    #print txt
    #rex6 = re.compile('doi:\s+(?P<doi>.+)').search(txt)
    #if rex: doi = rex.group('doi')
    #if rex2: doi = rex2.group('doi')
    #if rex3: doi = rex3.group('doi')
    #if rex4: doi = rex4.group('doi')
    #if rex5: doi = rex5.group('doi')
    #if rex6: doi = rex6.group('doi')
    doi = ''
    doi_line_re = re.compile('(?P<doi>.*DOI.*)').search(text)
    if doi_line_re:
        doi_line = doi_line_re.group('doi')
        if verbose:
            print 'DOI:', doi_line
        rrr = 'DOI:{0,1}\s{0,1}(?P<doi>[\w\d\.\-\\\/\–]+)'
        rex = re.compile(rrr).search(doi_line)
        doi = rex.group('doi')
        if verbose:
            print 'doi - found: ', doi
    return doi


def get_value(field, txt, verbose=False):
    """ """
    c = '(?P<line><meta .+"' + field + '".+>)'
    rex = re.compile(c).search(txt)
    if rex:
        line = rex.group('line')
        if verbose:
            print line
        b = 'content="(?P<pmid>\d+)"'
        bb = re.compile(b).search(line)
        return bb.group('pmid')


def get_pmid_via_doi_net(doi):
    """Give PMID using give DOI."""
    # here you must read in a field from the form
    params = urllib.urlencode({'hdl': doi})
    # the link where to send the data
    f = urllib.urlopen("http://dx.doi.org/", params)
    content = f.read()
    return get_value('citation_pmid', content)


def get_title_via_doi(doi, reference, customed_title, verbose=1):
    """Search in pubmex for a paper based on give doi.
    Get the paper, fetch pmid, get_title_via_pmid.

    input:
    * doi: e.g.
      and so on..
    * v(erbose)

    output:
    * title: string
    """
    doi = doi.replace('DOI:', '')
    doi = doi.replace('doi:', '')
    if DEBUG:
        print doi
    pmid = doi2pmid(doi)
    if DEBUG:
        print 'PMID:', pmid
    if pmid:
        return get_title_via_pmid(pmid, reference, customed_title)
    else:
        print 'ERROR: \t\tNot found in PubMed, although DOI (' + doi + ') was detected in the pdf!'
        pmid = get_pmid_via_doi_net(doi)
        return get_title_via_pmid(pmid, reference, customed_title)


def get_options(verbose=False):
    """ display options """
    usage = "%prog [options] id"
    description = """
examples: pubmex.py -p 17123955; pumex.py
-p 10.1038/embor.2008.212; pubmex.py -a -f file.pdf -r
"""
    version = __version__

    parser = optparse.OptionParser(description=description,
                              version=version,
                              usage=usage)

    parser.add_option("--PMID_or_DOI", "-p",
                      help="pass PMID/DOI of the paper")

    parser.add_option("--filename", "-f",
                      help="filename of a pdf")

    parser.add_option("--automatic", "-a", action="store_true",
                      help="try to get DOI automatically from a pdf, this option DOES NOT RENAME, use -r to force renaming")

    parser.add_option("--keywords", "-k", type="string",
                      help="""pass your keywords which makes a filename in a format 'RNA, structure' """)

    # TODO
    #parser.add_option("--reference", "-r", action="store_true",
    #                  help = "reference format", default=False)

    parser.add_option("--rename", "-r", action="store_true",
                      help="DOES rename files (only in a automatic mode)", default=False)

    if verbose:
        #(<Values at 0xb7b4b4cc: {'pmid': '1212'}>, [])
        print parser.parse_args()

    (options, arguments) = parser.parse_args()
    if DEBUG:
        print options

    if not options.PMID_or_DOI and not options.automatic:
        parser.print_help()
        parser.error('You have to pass PMID or DOI. You migth set -a (automatic)')
        sys.exit(1)

    return options, arguments


def rename(src, dst, rename_flag):
    if rename_flag:
        print 'RENAME:\t\t', src, '-*DO*->', dst
        shutil.move(src, dst)
    else:
        print '(fake, use -r)\t', src, '-*NOT*->', dst
        pass


def test_if_pdftotext():
    """Test if pdftotext is installed - quick-and-dirty hack."""
    return_code = subprocess.call("pdftotext -h", shell=True, stderr=subprocess.PIPE)
    if return_code == 99:
        pass
    else:
        print "Please, install pdftotext! Pdftotext is a part of the proper-utils package.\nFor debian-based linux systems run 'sudo apt-get install poppler-utils' and start pubmex again. Good luck!"
        sys.exit(1)

def main():
    test_if_pdftotext()
    OPTIONS, ARGUMENTS = get_options()
    keywords = ''
    title = ''
    # 1st mode: non-automatic
    if OPTIONS.PMID_or_DOI:
        if is_it_pmid(OPTIONS.PMID_or_DOI):
            title = get_title_via_pmid(OPTIONS.PMID_or_DOI, False, OPTIONS.keywords)
        else:
            title = get_title_via_doi(OPTIONS.PMID_or_DOI, False, OPTIONS.keywords)  # OPTIONS.reference
        #if title:
        #    print title
        #else:
        #    print 'ERROR: \t\tProblem! Check your PMID/DOI!'
        print title
    # 2nd mode: automatic
    if OPTIONS.automatic:
        filename = OPTIONS.filename
        print 'FILENAME:\t', filename
        text = pdf2text(filename)
        title = get_title_auto_from_text(text, False, OPTIONS.keywords)
        if title:
            print 'TITLE: \t\t', title
            ## basename = os.path.basename(filename)
            dirname = os.path.dirname(filename)
            if dirname == '':
                dirname = '.' + dirname  # .//file if dirname equals ''
            src = filename
            if DONE_FOLDER_NAME:
                try:
                    os.mkdir(dirname + os.sep + DONE_FOLDER_NAME)
                except OSError:
                    pass
                dst = dirname + os.sep + DONE_FOLDER_NAME + os.sep + title
            else:
                dst = dirname + os.sep + title
            rename(src, dst, OPTIONS.rename)
        else:
            get_pmid_via_search_in_pubmex_line_by_line(text)


        #if title:
        #    print title
        #else:
        #    print 'ERROR: \t\tProblem! The pubmex could not find DOI in the pdf file!'
        # TODO



if '__main__' == __name__:
    main()
