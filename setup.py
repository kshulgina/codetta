from setuptools import setup

setup(
    name="codetta",
    version="1.1.0",
    description="A Python program for predicting the genetic code (codon table) from nucleotide sequence data",
    author="Yekaterina Shulgina",
    author_email="shulgina@g.harvard.edu",
    url="https://github.com/kshulgina/codetta",
    python_requires='>=3.5',
    install_requires=['numpy>=1.18', 'scipy>=1.4'],
    py_modules = ['codetta', 'codetta_download', 'codetta_align', 'codetta_summary', 'codetta_infer', 'helper_functions']
)
