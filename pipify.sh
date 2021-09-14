#!/bin/zsh
trash dist
trash build
python setup.py bdist_wheel
twine upload dist/* --verbose
