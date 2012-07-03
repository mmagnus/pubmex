echo '# get a title based on PMID'
echo 'python pubmex.py -p 20533529 '
python pubmex.py -p 20533529 
echo
echo '# get data based on PMID, change a title to keywords'
echo \'python pubmex.py -p 20533529 -k 'MLH1,PMS2'\'
python pubmex.py -p 20533529 -k 'MLH1,PMS2'
echo
echo '# get data based on PMID, change a title to keywords + debug info'
echo \'python pubmex.py -d -p 20533529 -k 'MLH1,PMS2'\'
python pubmex.py -d -p 20533529 -k 'MLH1,PMS2'
echo
echo '# get data based on DOI'
echo 'python pubmex.py -p 10.1038/embor.2008.212'
python pubmex.py -p 10.1038/embor.2008.212
echo
echo '# get data based on a given pdf file, show a new new, DO NOT rename'
echo 'python pubmex.py -a -f demo/demo01.pdf'
python pubmex.py -a -f demo/demo01.pdf
echo
echo '# get data based on a given pdf file, show a new new, DO rename'
echo 'pubmex.py -a -f demo/demo01.pdf -r'
python pubmex.py -a -f demo/demo01.pdf -r

