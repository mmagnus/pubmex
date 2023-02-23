#! /usr/bin/env python3
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

import sys
import argparse
import subprocess

import os
import re
import tempfile
import subprocess
import shutil
import urllib.request, urllib.parse, urllib.error

from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='', includeContext=True)


MAIL = 'your_mail@gmail.com'
JDICT = { 'NUCLEIC.ACIDS.RES': 'NAR', 'NucleicAcidsRes' : 'NAR', 'BiochimBiophysActa': 'BBA',     }
# FirstAuthor.LastAuthor.keywords.Journal.year.pdf
DELIMITER_AUTHOR_TITLE = '-'
DELIMITER_KEYWORDS_JOURNAL = '.'
ADD_PMID = False
WORDS_TO_REMOVE = 'a, as, at, for, from, he, her, his, if, in, it, its, of, on, she, so, the, their, them, they, to, which' + ',with, and, by, during'
DONE_MOVE_TO_FOLDER = False
DONE_FOLDER_NAME = 'done'
HOW_MANY_LINES_TO_READ = 10
LENGHT_OF_LINE = 20
TEMPFILE_NAME = 'temp'
LJUST = 20
LJUST_SPACER = '.'
###############################
###############################


def hr():
    """Draw -------------."""
    ic('-' * 80)

def exe(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err

def dot(text):
    """Change ' '(space) into . ,e.g. Marcin Magnus -> Marcin.Magnus."""
    return text.replace(' ', '.')


def clean_string(text):
    """Execute operations to clean up a given string."""
    for word in WORDS_TO_REMOVE.split(','):
        text = text.replace('.' + word.strip() + '.', '.')
        text = text.replace('.' + word.strip().upper() + '.', '.')
    text = text.strip()
    text = text.replace(' ', '.')
    text = text.replace('"', '')
    text = text.replace(':', '.')    
    text = text.replace('?', '.')
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


# !!experimetal!! 
def get_pmid_via_search_in_pubmex_line_by_line(
        text='RNA tertiary structure prediction with ModeRNA.', debug=False):
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
    #    print txtfn, "is going to be opend"
    c = 0
    return None
    if not text:
        return None
    lines = text.split('\n')
    for line in lines:
        if (line.strip()) and (len(line) > LENGHT_OF_LINE):
            query = line
            if debug:
                ic(query)
            ic(query)
            ###
            if 0:
                from Bio import Entrez
                result = Entrez.esearch(db="pubmed", term=query, email=MAIL)
                d = Entrez.read(result)
                #print 'd[Count]: ', d['Count']
                pmid = d['IdList']
                if pmid:
                    if debug:
                        ic('pmid: ', pmid)
                if len(pmid) == 1:  ## uniq PMID is found
                    return pmid
                else:
                    for p in pmid:
                        if debug:
                            ic(p, get_title_via_pmid(p, debug))
                if c > HOW_MANY_LINES_TO_READ:
                    break
                c += 1

    y = 0
    for y in range(0, HOW_MANY_LINES_TO_READ):
        query = lines[y] + ' ' + lines[y+1]
        if debug:
                ic('query: ', query)
        result = Entrez.esearch(db="pubmed", term=query, email=MAIL)
        d = Entrez.read(result)
        #print 'd[Count]: ', d['Count']
        pmid = d['IdList']
        if pmid:
            ic('pmid: ', pmid)
            if len(pmid) == 1:  ## uniq PMID is found
                return pmid
    

def get_title_via_pmid(pmid, debug, reference='', customed_title=''):
    """Use biopython to get summary dict.

    input:
     * pmid of a publication, e.g. 17123955
       and so on..
     * v(erbose)

    output:
     * formated title(string)

  <title>
   Structure and mechanism of a methyltransferase ribozyme
  </title>
  
    """
    import urllib.request
    from bs4 import BeautifulSoup
    if pmid:
        ic(pmid)
        o = urllib.request.urlopen('http://pubmed.ncbi.nlm.nih.gov/' + pmid)
        txt = o.read()
        soup = BeautifulSoup(txt, 'html.parser')

        # if debug: ic(soup.prettify())
    
        title = soup.title.string
        # content="Nat Chem Biol" name="citation_publisher">
        journal = re.findall(r'content="(.+)" name="citation_publisher"', soup.prettify())[0]
        ic(journal)
        journal = journal.replace(' ', '')  # Nat Chem Biol -> NatChemBiol
        # <meta content="2022/03/17" name="citation_online_date">
        year = re.findall(r'meta content="(\d+)/\d+/\d+" name="citation_online_date', soup.prettify())
        if not year:
            # <meta content="2007/3" name="citation_publication_date">
            #  <meta content="2016/11/29" name="citation_publication_date">
            year = re.findall(r'meta content="(\d+)/.+" name="citation_publication_date', soup.prettify())
        year = year[0]
        ic(year)
        """
> auths: ['Jie Deng',
          'Timothy J Wilson',
          'Jia Wang',
          'Xuemei Peng',
          'Mengxiao Li',
          'Xiaowei Lin',
          'Wenjian Liao',
          'David M J Lilley',
          'Lin Huang']
        """
        auths = re.findall(r'meta content="(.+)" name="citation_author"', soup.prettify())# [0]
        ic(auths)
        
        first = auths[0].split()[-1]
        last = auths[-1].split()[-1] # take last name
        ic(first, last)
        if first == last:
            last = ''
        ic(first, last)
        all = first + '.' + last + '-' + title + '-' + journal + '.' + year + '.pdf'
        all = all.replace(' ', '.')
        return all

    return None

    result = Entrez.esearch(db="pubmed", term=pmid, email=MAIL) # id
    ic(Entrez.read(result))
    try:
        summary_dict = Entrez.read(result)[0]
    except RuntimeError:  # there is NO pmd like this
        return False

    if debug:
        ic('summary_dict'.ljust(LJUST, LJUST_SPACER), summary_dict)

    if not reference:
        title_of_pub = dot(summary_dict['Title'].strip())

        first_author = summary_dict['AuthorList'][0].split(' ')[0].strip()
        title = first_author

        last_author = summary_dict['LastAuthor'].split(' ')[0].strip()
        if first_author != last_author:
            title += '.' + last_author

        if not customed_title:
            title += DELIMITER_AUTHOR_TITLE + title_of_pub
        else:
            title += DELIMITER_AUTHOR_TITLE + prepare_customed_title(customed_title)

        if ADD_PMID:
            title += DELIMITER_AUTHOR_TITLE + pmid

        title += DELIMITER_AUTHOR_TITLE + dot(summary_dict['Source'].replace(' ', ''))
        title += '-' + str(summary_dict['PubDate'].split()[0])
        title += '.pdf'

        title = clean_string(title)

        for k in list(JDICT.keys()):
            title = title.replace(k, JDICT[k])

    # TODO refernece
    else:
        title = ", ".join(summary_dict['AuthorList']) + " " + summary_dict['Title'].strip() + " " + summary_dict['Source'].strip() + " " + summary_dict['SO']
    return title


def text2clip(title):
    """Output is sent to the clipboard."""
    cmd = "echo '" + title.strip() + "'| xclip -selection clipboard"
    os.system(cmd)


def doi2pmid(doi, debug):
    """Get a PMID based on a given DOI.

    https://www.ncbi.nlm.nih.gov/esearch.fcgi?db=pubmed&term=asthma&field=title
    """
    # *needto be improved*
    # doi = doi.replace('/', '-')
    import urllib.request
    o = urllib.request.urlopen('http://pubmed.ncbi.nlm.nih.gov/?term=%20' + doi + '&sort=date')
    pmid = re.compile('data-article-pmid="(?P<pmid>\d+)"').search(str(o.read()))#soup.prettify())
    if pmid:
        return pmid['pmid']
    return None

def title2pmid(title, debug):
    """Get a PMID based on a given DOI.

    https://www.ncbi.nlm.nih.gov/esearch.fcgi?db=pubmed&term=asthma&field=title
    """
    # *needto be improved*
    # doi = doi.replace('/', '-')
    import urllib.request
    from urllib.parse import quote
    url = 'http://pubmed.ncbi.nlm.nih.gov/?term=' + quote(title) + '&sort=date'
    ic(url)
    o = urllib.request.urlopen(url)
    pmid = re.compile('data-article-pmid="(?P<pmid>\d+)"').search(str(o.read()))#soup.prettify())
    ic('xx', pmid)
    if pmid:
        return pmid['pmid']
    return None


def get_title_auto_from_text(text, debug, reference, customed_title):
    doi = get_doi_from_text(text, debug)
    if doi:
        title = get_title_via_doi(doi, debug, reference, customed_title)
    if not doi or not title:
        if debug: ic('DOI has *not* been found automatically!')
        pmid = get_pmid_via_search_in_pubmex_line_by_line(text, debug)
        title = get_title_via_pmid(pmid, debug)
    return title
    
def pdf2text(filename, debug):
    """Convert a pdf file to a flat file.

    input:
    * filename

    caution:
    * pdftotext: is required
    """
    global TEMPFILE_NAME
    if debug:
        ic(('generate ./' + TEMPFILE_NAME).ljust(LJUST, LJUST_SPACER) + '[OK]')
    else:
        f = tempfile.NamedTemporaryFile()
        TEMPFILE_NAME = f.name
    txtfn = TEMPFILE_NAME
    cmd = 'pdftotext ' + filename + ' ' + txtfn

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    out = process.stdout.read().decode()
    err = process.stderr.read().decode()

    if debug:
        ic('out:', out)
        ic('err:', err)

    if not err:
        txt = open(txtfn, encoding="utf8", errors='ignore').read()
        if debug:
            ic(txtfn, "is going to be opened")
        return txt
    else:
        ic('ERROR: pdftotext ' + filename + ' ' + txtfn)
    
def get_doi_from_text(text, debug):
    """Use several regular expression are used to get DOI in a text.

    input:
    * text
    * debug
    
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

    if not text:
        return None
    text = text.upper()
    #text = text.replace('–', '-')
    #print txt
    #rex6 = re.compile('doi:\s+(?P<doi>.+)').search(txt)
    #if rex: doi = rex.group('doi')
    #if rex2: doi = rex2.group('doi')
    #if rex3: doi = rex3.group('doi')
    #if rex4: doi = rex4.group('doi')
    #if rex5: doi = rex5.group('doi')
    #if rex6: doi = rex6.group('doi')
    doi_line_re = re.compile('(?P<doi>.*DOI.*)').search(text)
    if doi_line_re:
        doi_line = doi_line_re.group('doi')
        if debug:
            ic('doi_line: '.ljust(LJUST, LJUST_SPACER), doi_line)
        # fixes
        doi_line = doi_line.replace('10. 1073', '10.1073') # for pnas

        rrr = 'DOI:{0,1}\s{0,1}(?P<doi>[\w\d\.\-\\\/\–]+)'
        rex = re.compile(rrr).search(doi_line)
        if rex:
            doi = rex.group('doi')
            # .ORG/
            if doi.startswith('.ORG'):
                doi = doi.replace('.ORG/','')
        if debug:
            ic('doi is found: '.ljust(LJUST, LJUST_SPACER), doi)
        return doi
    return None

def get_value(field, txt, debug):
    """ """
    c = '(?P<line><meta .+"' + field + '".+>)'
    rex = re.compile(c).search(txt)
    if rex:
        line = rex.group('line')
        if debug:
            ic(line)
        b = 'content="(?P<pmid>\d+)"'
        bb = re.compile(b).search(line)
        return bb.group('pmid')


def get_pmid_via_doi_net(doi):
    """Give PMID using give DOI."""
    # here you must read in a field from the form
    params = urllib.parse.urlencode({'hdl': doi})
    # the link where to send the data
    try:
        f = urllib.request.urlopen("http://dx.doi.org/", params)
        content = f.read()
    except:
        return ''
    else:
        return get_value('citation_pmid', content)


def get_title_via_title(title, debug):
    return ''
    pmid = title2pmid(title, debug)
    ic(pmid)
    return get_title_via_pmid(pmid, debug)#, reference, customed_title)

def get_title_via_doi(doi, debug=False, reference='', customed_title=''):
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
    if len(doi) < 10:
        return ''
        
    if debug:
        ic(doi)
    pmid = doi2pmid(doi, debug)
    if debug:
        ic(pmid)
    if pmid:
        return get_title_via_pmid(pmid, debug, reference, customed_title)
    else:
        if debug:
            ic('Not found in PubMed, although DOI (' + doi + ') was detected in the pdf!')
        pmid = get_pmid_via_doi_net(doi)
        return get_title_via_pmid(pmid, debug, reference, customed_title)

def get_parser():
    """ display options """
    usage = "%prog [options] id"
    description = """
examples: pubmex.py -p 17123955; pumex.py
-p 10.1038/embor.2008.212; pubmex.py -a -f file.pdf -r
"""

    parser = argparse.ArgumentParser(
        description='', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")

    parser.add_argument("file", help="", default="", nargs='+')

    parser.add_argument("-p", "--PMID_or_DOI",
                      help="pass PMID/DOI of the paper")

    #parser.add_argument("-f", "--file",
    #                    help="filename (files) of a pdf", nargs='+')

    parser.add_argument("-k", "--keywords", 
                      help="""pass your keywords which makes a filename in a format 'RNA, structure' """)

    ## # TODO
    ## #parser.add_option("--reference", "-r", action="store_true",
    ## #                  help = "reference format", default=False)

    parser.add_argument("-d", "--debug", action="store_true",
                       help="show debug message", default=False)

    return parser


def rename(src, dst, rename_flag):
    # if debug: ic(rename_flag)
    if rename_flag:
        ic('mv ', src, '-->', dst)
        shutil.move(src, dst)
    else:
        ic('CAUTION! THE FILE WAS NOT RENAME\n', 'mv', src, '-->', dst)
        pass


def test_if_pdftotext():
    """Test if pdftotext is installed - quick-and-dirty hack."""
    return_code = subprocess.call("pdftotext -h", shell=True, stderr=subprocess.PIPE)
    if return_code == 0 or return_code == 99:
        pass
    else:
        ic("Please, install pdftotext! Pdftotext is a part of the proper-utils package.\nFor debian-based linux systems run 'sudo apt-get install poppler-utils' and start pubmex again. Good luck!")
        sys.exit(1)


# main
def main():
    test_if_pdftotext()

    parser = get_parser()
    args = parser.parse_args()

    title = ''
    debug = args.debug
    if not debug:
        ic.disable()

    # 1st mode: non-automatic
    if args.PMID_or_DOI:
        if is_it_pmid(args.PMID_or_DOI):
            title = get_title_via_pmid(args.PMID_or_DOI, debug, False, args.keywords)
        else:
            title = get_title_via_doi(args.PMID_or_DOI, debug, False, args.keywords)  # opt.reference
        #if title:
        #    print title
        #else:
        #    print 'ERROR: \t\tProblem! Check your PMID/DOI!'
        if debug: ic(title)

    # 2nd mode: automatic
    if args.file:
        #################################
        if list != type(args.file):
            args.file = [args.file]
        ##################################
        for filename in args.file:
            # if osx
            ic('-' * 80)
            
            ic(filename)

            # import pdf2doi
            import pdf2doi
            pdf2doi.config.set('verbose', False)
            doi =  pdf2doi.pdf2doi(filename)
            title = get_title_via_doi(doi['identifier'], debug=debug, reference='', customed_title='')
            # ic.disable()
            if not title:
                # check if filename is a DOI
                fn = os.path.basename(filename)    
                title = get_title_via_doi(fn.replace('.pdf', '').replace(':', '/'), args.debug, False, args.keywords)
                title_from_tag = ''

            if not title:
                from sys import platform
                if platform == "darwin":
                    out, err = exe('mdls ' + filename)
                    for l in out.split('\n'):
                        if 'kMDItemTitle' in l:
                            ic(l)
                        # kMDItemTitle
                        # kMDItemDescription                     = "Nature Chemical Biology, doi:10.1038/s41589-022-00982-z"
                        if 'doi' in l:  # kMDItemDescription      = "Nature Cell Biology 18, 1261 (2016). doi:10.1038/ncb3446"
                            r = re.search('doi:(?P<doi>.+)', l)
                            if r:
                                doi = r.group('doi').replace('"','')
                                title = get_title_via_doi(doi, args.debug, False, args.keywords)
                                ic(doi, title)


                    if not title:  # if not yet title
                        for l in out.split('\n'):
                            if 'kMDItemTitle' in l:
                            # kMDItemTitle                           = "CSSR: assig"
                                tag, title_from_tag = l.split('=')
                                title_from_tag = clean_string(title_from_tag) + '.pdf' # .replace(' ', '.')
                                # title = get_title_via_title(title, args.debug)#, False, args.keywords)
                                # ic(title)
                                
            if 1:  # search for doi in text
                if not title:  # if title not yet found
                    text = pdf2text(filename, args.debug)
                    title = get_title_auto_from_text(text, args.debug, False, args.keywords)

            if not title:
                if title_from_tag:
                    title = title_from_tag
                
            if title:
                title = clean_string(title)
                if debug:
                    ic('the title is ... ', title)
                dirname = os.path.dirname(filename)
                if dirname == '':
                    dirname = '.' + dirname  # .//file if dirname equals ''
                src = filename
                #if DONE_FOLDER_NAME:
                #    try:
                #        os.mkdir(dirname + os.sep + DONE_FOLDER_NAME)
                #    except OSError:
                #        pass
                #    dst = dirname + os.sep + DONE_FOLDER_NAME + os.sep + title
                #else:
                dst = dirname + os.sep + title
                rename(src, dst, not args.debug)  # rename)
            else:
                ic('ERROR: \t\tProblem! The pubmex could not find automatically a title for the pdf file! Sorry!')

            print(title)
            
if '__main__' == __name__:
    main()
