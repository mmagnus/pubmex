# Pu(b)mex

Marcin Magnus (m.magnus@o2.pl)

# DESCRIPTION

`pubmex.py` is a script to get a fancy paper title based on given DOI or PMID.

Format of the title:

	a first author . a last author - (title("dotted") or your customed title) . PMID . journal . year . pdf

Using xclip pubmex.py automatically copy the title to the clipboard. Just type Ctrl+V or paste it somewhere else.

# DEPENDENCIES

- biopython (http://biopython.org/wiki/Biopython)
- pdftotext (http://linux.die.net/man/1/pdftotext)
- xclip (http://sourceforge.net/projects/xclip/) (not required)

# INSTALLATION

## Ubuntu (Debian-based system)

	apt-get install xclip python-biopython pdftotext

.. clone & run it :-) Set PATHs to `pubmex.py`

## MAC OSX

    sudo port install poppler #pdftotext
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

Visit https://github.com/m4rx9/pubmex/wiki

Run a script with examples:

    ./examples-test.sh

# BUGS

`pubmex.py` might have still problem if it has to parse title of authors with characters like ółą and so on.

# COPYRIGHT AND LICENCE

pubmex.py is Copyright (C) 2010-2015 Marcin Magnus.  All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as GLP

# AUTHOR INFORMATION

Marcin Magnus, m.magnus@o2.pl

# ACKNOWLEDGEMENTS

@todo

# HISTORY

version 0.3

- osx installation

version 0.2

- small changes




