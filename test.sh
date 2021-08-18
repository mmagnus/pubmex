# echo '# get a title based on PMID'
# echo 'python pubmex.py -p 20533529 '
# python pubmex.py -p 20533529 
# echo '# get data based on PMID, change a title to keywords'
# echo \'python pubmex.py -p 20533529 -k 'MLH1,PMS2'\'
# python pubmex.py -p 20533529 -k 'MLH1,PMS2'
# echo '# get data based on PMID, change a title to keywords + debug info'
# echo \'python pubmex.py -d -p 20533529 -k 'MLH1,PMS2'\'
# python pubmex.py -d -p 20533529 -k 'MLH1,PMS2'
# echo '# get data based on DOI'
# echo 'python pubmex.py -p 10.1038/embor.2008.212'
# python pubmex.py -p 10.1038/embor.2008.212
set -x
python pubmex.py  demo/demo01.pdf --debug

python pubmex.py demo/demo01.pdf --debug

python pubmex.py demo/ncb3446.pdf --debug

python pubmex.py demo/10.1261:rna.418407.pdf --debug
