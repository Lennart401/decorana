# DECORANA Python Wrapper
 
This library provides Python wrappers for the classic program for ecological data analyses, DECORANA (_DEtrended COrrespondence ANAlysis_). The program was written by M. O. Hill in FORTRAN for mainframe computers, and modified for the IBM PC. 

This modified version uses the "strict" convergence criteria of Oksanen & Minchin (1997) for eigenanalysis, with a tolerance of 0.000005 and a maximum iteration limit of 999. In DECORANA, the bug in non-linear scaling has been corrected.

This library ships with the original fortran source code of the DECORANA program, which is compiled every time the decorana function is called. Thus, the GNU Fortran Compiler (gfortran) is required to be installed on your system and available in the path. The wrapper checks for the availability of gfortran and raises an error if it is not found. You can get the fortran compiler from [here](https://gcc.gnu.org/wiki/GFortranBinaries). For ubuntu, you can install it as package using the command `sudo apt install gfortran`. 

The fortran program requires in input data file in Cornell Condensed format, containing the community data to be analysed. The layout of this file should follow the same rules as for the original versions of DECORANA. The __write_cep() function included in this package will automatically convert a pandas dataframe into the required format for input to the program. More information can be found here: http://ordination.okstate.edu/formats.htm.

Most of this library has been adopted from the CornPy package (https://github.com/maurobio/cornpy), and the biplot function has been adopted from the ecopy package (https://github.com/Auerilas/ecopy). The csv files under the /data directory have been copied from the CornPy package.

# Installation

1. Install the GNU Fortran Compiler (gfortran) from [here](https://gcc.gnu.org/wiki/GFortranBinaries). Make sure the gfortran binary is available in the PATH.
2. Clone this repository

    `git clone https://github.com/Lennart401/decorana.git`
3. Install the package using pip

    `pip install .`

# License
**DECORANA** is distributed under the GNU General Public License

# Version
0.1.0

# Examples

    import pandas as pd
    from decorana import decorana, biplot

    df = pd.read_csv('./data/gauch.csv', index_col=0)
    site_scores, species_scores, site_labels, species_labels = decorana(df)
    biplot(site_scores, species_scores, site_labels, species_labels)

