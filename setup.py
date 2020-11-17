import setuptools
import glob

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name="airpg",
	version="0.1.5",
	author="Tilman Mehl, Michael Gruenstaeudl",
	author_email="tilmanmehl@zedat.fu-berlin.de, m.gruenstaeudl@fu-berlin.de",
    description="A package to automatically access the inverted repeats of archived plastid genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/michaelgruenstaeudl/airpg',
    download_url='https://github.com/michaelgruenstaeudl/airpg/archive/v0.1.5.tar.gz',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: BSD License',
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
    python_requires='>=3.6',
    keywords='plastid genomes, inverted repeats, NCBI Nucleotide',
    license='BSD',
    packages=['airpg'], # So that the subfolder 'airpg' is read immediately.
    #packages = find_packages(),
    install_requires=['biopython', 'ete3', 'entrezpy', 'pandas'],
    scripts=['airpg/scripts/airpg_retrieve.py', 'airpg/scripts/airpg_analyze.py', 'airpg/scripts/airpg_generate_blacklist.py'],
    test_suite='setup.my_test_suite',
    include_package_data=True,
    zip_safe=False
)
