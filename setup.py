import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name="airpg",
	version="0.1.8",
	author="Tilman Mehl, Michael Gruenstaeudl",
	author_email="tilmanmehl@zedat.fu-berlin.de, m.gruenstaeudl@fu-berlin.de",
    description="A package to automatically access the inverted repeats of archived plastid genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/michaelgruenstaeudl/airpg',
    download_url='https://github.com/michaelgruenstaeudl/airpg/archive/v0.1.8.tar.gz',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
    python_requires='>=3.6',
    keywords='plastid genomes, inverted repeats, NCBI Nucleotide',
    license='GPLv3',
    #packages=['airpg'], # So that the subfolder 'airpg' is read immediately.
    packages = setuptools.find_packages(),
    #entry_points={
    #  "console_scripts": ["airpg-retrieve=airpg.scripts.airpg_retrieve", "airpg-analyze=aripg.scripts.airpg_analyze", "airpg-make-blacklist=airpg.scripts.airpg_generate_blacklist"]
    #},
    install_requires=['biopython', 'ete3', 'entrezpy', 'pandas'],
    scripts=['airpg/scripts/airpg_retrieve.py', 'airpg/scripts/airpg_analyze.py', 'airpg/scripts/airpg_generate_blacklist.py'],
    test_suite='setup.my_test_suite',
    include_package_data=True,
    zip_safe=False
)
