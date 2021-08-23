# Pu(b)mex

[![tag](https://img.shields.io/github/release/mmagnus/pubmex.svg)](https://github.com/mmagnus/pubmex/releases)

`pubmex.py` is a script to get a fancy paper title based on given DOI or PMID.

Format of the title:

	a first author . a last author - (title("dotted") or your customed title) . PMID . journal . year . pdf
	e.g.
	  Kelley.Scott.The.evolution.biology.shift.towards.engineering.prediction-generating.tools.away.traditional.research.practice.EMBORep.2008.pdf


Nowadays, it’s not a big issue, with all Mendeley and other tools, however...

I don’t want to put any PDF file collected on the way into my library, because then it gets super big (and then it’s hard to sync it for example with Dropbox). So now I can keep these PDF files into pdf-icebox and re-name them niecely automatically:

```
$ ls
Hnisz.Sharp.Phase.Separation.Model.Transcriptional.Control.Cell.2017.pdf
Sharp.Hockfield.Convergence.The.future.health.Science.2017.pdf
```

Usage:

```
    $ pubmex.py sharp2017.pdf
    Sharp.Hockfield.Convergence.The.future.health.Science.2017.pdf
    mv  sharp2017.pdf --> ./Sharp.Hockfield.Convergence.The.future.health.Science.2017.pdf

    $ pubmex.py Query.Konarska.pdf
    the title is ...  Smith.Konarska."Nought.may.endure.but.mutability".spliceosome.dynamics.regulation.splicing.MolCell.2008.pdf
    
    $ pubmex.py eabc9191.full.pdf
    mv  eabc9191.full.pdf --> ./Balas.Johnson.Establishing.RNA-RNA.interactions.remodels.lncRNA.structure.promotes.PRC2.activity.SciAdv.2021.pdf
```

# DEPENDENCIES

- biopython (http://biopython.org/wiki/Biopython)
- pdftotext (http://linux.die.net/man/1/pdftotext)

# INSTALLATION

Python 3 is required.

    git clone https://github.com/mmagnus/pubmex
    pip install -e pubmex

## Ubuntu (Debian-based system)

	apt-get install xclip python-biopython pdftotext

.. clone & run it :-) Set PATHs to `pubmex.py`

## MAC OSX

    sudo port install poppler #pdftotext
	sudo port install py-pip # 
	sudo pip install biopython # or sudo port install biopython # (but it did't work for me)
	
.. clone & run it :-) Set PATHs to `pubmex.py`

## WINDOWS (NOT TESTED!!!!!!!)

Install:

- Python http://www.python.org/ftp/python/2.6.4/python-2.6.4.msi
- Biopython http://biopython.org/DIST/biopython-1.53.win32-py2.6.exe
- GTK+ http://sourceforge.net/projects/gtk-win/ direct link http://downloads.sourceforge.net/project/gtk-win/GTK%2B%20Runtime%20Environment/GTK%2B%202.16/gtk2-runtime-2.16.6-2010-02-24-ash.exe?use_mirror=garr
- pyGTK http://ftp.gnome.org/pub/GNOME/binaries/win32/pygtk/2.16/pygtk-2.16.0+glade.win32-py2.6.exe
- and pdftotext if you want to use `-a`

.. and it might work. I don't know (I don't care ;-))

# EXAMPLES

Visit https://github.com/mmagnus/pubmex/wiki

Run a script with examples:

    ./examples-test.sh

# BUGS

`pubmex.py` might have still problem if it has to parse title of authors with characters like ółą and so on.

# HISTORY

version 1.1 simplify input, pubmex.py *.pdf

version 1.0: with recent bugfixes 2021

version 0.3

- osx installation

version 0.2

- small changes
