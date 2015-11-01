from setuptools import setup, find_packages


def get_version():
    f = open('pubmex.py')
    try:
        for line in f:
            if line.startswith('__version__'):
                return eval(line.split('=')[-1])
    finally:
        f.close()


def get_long_description():
    descr = []
    for fname in ['README.md']:  # , # 'CHANGES.txt':   # , 'TODO.txt'
        f = open(fname)
        try:
            descr.append(f.read())
        finally:
            f.close()
    return '\n\n'.join(descr)


setup(
    name='pubmex',
    version=get_version(),
    description="pu(b)mex: a scientific publication renamer",
    long_description=get_long_description(),
    keywords='publication, pubmed, doi, pmid',
    author='Marcin Magnus',
    author_email='m.magnus@o2.pl',
    url='https://github.com/m4rx9/pubmex/wiki',
    license='GPLv3',
    py_modules=['pubmex'],
    namespace_packages=[],
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        #'setuptools',
        # -*- Extra requirements: -*-
    ],
    data_files=[('', ['README.md']),
                ('demo', ['demo/demo01.pdf', 'demo/demo02.pdf'])
                ],
    entry_points={
        'console_scripts': [
            'pubmex = pubmex:main',
        ],
    },
)
