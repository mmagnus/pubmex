from setuptools import setup, find_packages
import versioneer

# read the contents of your README file
with open('README.md') as f:
    long_description = f.read()
    
setup(
    name='pubmex',
    description="pu(b)mex: a scientific publication renamer",
    long_description=long_description,
    version='1.4.2', # versioneer.get_version(),
    # cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    test_suite="tests",
    url='https://github.com/mmagnus/pubmex',
    scripts=['pubmex.py',
             ],
    license='MIT',
    author='Marcin Magnus',
    author_email='mag_dex@o2.pl',
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords='publication, pubmed, doi, pmid',
    install_requires=[
        'icecream'
       ],
)
