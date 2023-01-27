import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = [
    'numpy>=1.24.1',
    'matplotlib>=3.6.3',
    'pandas>=1.5.3',
]

setuptools.setup(
    name="decorana",
    version="0.1.0",
    author="Lennart Popkes",
    author_email="lennart.popkes@gmail.com",
    description="DECORANA in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Lennart401/decorana",
    install_requires=install_requires,
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=['ordination', 'classification', 'ecology', 'multivariate data analysis', 'DECORANA',
              'detrended correspondence analysis'],
    package_data={
        'decorana': ['data/*.csv'],
    },
)