import numpy as np
import pandas as pd
import shutil
import subprocess
import os
import matplotlib.pyplot as plt

from pathlib import Path


def decorana(data, iweight=0, iresc=0, ira=0, mk=0, short=0):
    """Detrended correspondence analysis and basic reciprocal averaging.

    This function performs detrended correspondence analysis and basic
    reciprocal averaging or orthogonal correspondence analysis. It is a
    Python wrapper around the program DECORANA written in Fortran by
    M.O. Hill (1979).

    Parts of this function has been adopted from https://github.com/maurobio/cornpy.

    Args:
        data (pd.DataFrame): A pandas dataframe
        iweigh (int): Downweighting of rare species (default = 0: no)
        iresc (int): Number of rescaling cycles (default = 0: no rescaling)
        ira	(int): Type of analysis (0: detrended, 1: basic reciprocal averaging, default = 0)
        mk (int): Number of segments in rescaling (default = 0)
        short (int): Shortest gradient to be rescaled (default = 0)

    Returns:
        site_scores (np.ndarray): array of site scores
        species_scores (np.ndarray): array of species scores
        site_labels (list): list of site labels
        species_labels (list): list of species labels
    """
    # compile program
    executable_path = __compile_fortran()

    # write params.dat and cep.dat
    nrows, ncols, row_lab, col_lab = __write_cep(data)
    with open('params.dat', 'w') as parfile:
        content = f"""cep.dat
        -1
        0
        {iweight}
        {iresc}
        {ira}
        {mk}
        {short}
        1
        0
        """
        parfile.write(content)

    # call decorana executable
    process_result = subprocess.run(f'{executable_path} < params.dat',
                                    shell=True,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL)
    if process_result.returncode != 0:
        raise Exception(f'decorana.exe exited with code {process_result.returncode}!')

    # remove decorana executable, params.dat and cep.dat
    os.remove('decorana.exe')
    os.remove('cep.dat')
    os.remove('params.dat')

    # read site and species scores
    site_scores = np.genfromtxt('decorana.out', usecols=(0, 1, 2, 3), max_rows=nrows)
    species_scores = np.genfromtxt([s[41:] for s in open('decorana.out')], usecols=(0, 1, 2, 3), max_rows=ncols)

    # remove decorana.out and decorana.prt
    os.remove('decorana.out')
    os.remove('decorana.prt')

    # return scores and labels
    return site_scores, species_scores, row_lab, col_lab


def biplot(siteScores, spScores, siteLabs=None, spLabs=None, xax=1, yax=2, showSp=True, showSite=True, spCol='r',
           siteCol='k', spSize=6, siteSize=6, markersize=4, xlim=None, ylim=None):
    """Biplot sites and species from DCA.

    This function has been adopted from ...

    Args:
        siteScores (np.ndarray): site scores
        spScores (np.ndarray): species scores
        siteLabs (list): labels of the sites, if None, the labels will not be plotted
        spLabs (list): labsl of the species, if None, the labels will not be plotted
        xax (int): x-Axis
        yax (int): y-Axis
        showSp (bool): True to show species
        showSite (bool): True to show sites
        spCol (str): color to plot species in
        siteCol (str): color to plot sites in
        spSize (int): size of species labels
        siteSize (int): size of site labels
        markersize (int): size of the markers
        xlim (int): limit of x axis
        ylim (int): limit of y axis

    Returns:
        None
    """
    siteScores = pd.DataFrame(siteScores)
    spScores = pd.DataFrame(spScores)

    f, ax = plt.subplots()
    if showSite:
        ax.plot(siteScores.iloc[:, xax - 1], siteScores.iloc[:, yax - 1], f'{siteCol}o',
                ms=markersize, label='Site Scores')
        if siteLabs is not None:
            [ax.text(x, y, s, fontsize=siteSize, color=siteCol, ha='center', va='bottom') for x, y, s in
             zip(siteScores.iloc[:, xax - 1], siteScores.iloc[:, yax - 1], siteLabs)]
    if showSp:
        ax.plot(spScores.iloc[:, xax - 1], spScores.iloc[:, yax - 1], f'{spCol}^',
                ms=markersize, label='Species Scores')
        if spLabs is not None:
            [ax.text(x, y, s, fontsize=spSize, color=spCol, ha='center', va='bottom') for x, y, s in
             zip(spScores.iloc[:, xax - 1], spScores.iloc[:, yax - 1], spLabs)]
        xmax = max(np.amax(siteScores.iloc[:, xax - 1]), np.amax(spScores.iloc[:, xax - 1]))
        xmin = min(np.amin(siteScores.iloc[:, xax - 1]), np.amin(spScores.iloc[:, xax - 1]))
        ymax = max(np.amax(siteScores.iloc[:, yax - 1]), np.amax(spScores.iloc[:, yax - 1]))
        ymin = min(np.min(siteScores.iloc[:, yax - 1]), np.min(spScores.iloc[:, yax - 1]))
        ax.set_xlim([xmin * 1.15, xmax * 1.15])
        ax.set_ylim([ymin * 1.15, ymax * 1.15])
    if xlim is not None:
        if not isinstance(xlim, list):
            msg = "xlim must be a list"
            raise ValueError(msg)
        ax.set_xlim(xlim)
    if ylim is not None:
        if not isinstance(ylim, list):
            msg = 'ylim must be a list'
            raise ValueError(msg)
        ax.set_ylim(ylim)
    ax.set_xlabel('CA Axis {!s}'.format(xax))
    ax.set_ylabel('CA Axis {!s}'.format(yax))
    ax.legend()
    plt.show()


def __compile_fortran(path='decorana.f') -> Path:
    """Compile fortran file to binary executable.

    Args:
        path (str): filename/path of the fortran source code file

    Return:
        (Path): path to the binary executable
    """
    # check if fortran compiler exists in path
    if (compiler := shutil.which('gfortran')) is not None:
        # compile script at path to decorana.exe
        cwd = Path(os.getcwd())
        executable_path = cwd / 'decorana.exe'
        process = subprocess.run([compiler, '-o', executable_path, path],
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL)
        if process.returncode == 0:
            return executable_path
        else:
            raise Exception(f'Compilation of {path} failed!')
    else:
        raise Exception('GNU Fortran compiler is not installed!')


def __write_cep(data):
    """Converts data matrix to the CEP compressed format.

    This function converts a data matrix from normal form, with species as rows
    and samples as columns, to the Cornell Ecology Programs condensed format.
    This condensed format consists of data ponts entered as couplets
    consisting of the number for the species and the abundance. Each line of
    the file begins with the number of the sample, followed by the couplets.
    The data for a sample may continue onto other lines. See the user's manuals
    for DECORANA and TWINSPAN for details on the structure of the CEP data files.

    This function has been adopted from https://github.com/maurobio/cornpy and
    adjusted to meet the Python3 standard and fixed a few issues which stopped
    the function from working correctly.

    Args:
        data (pd.DataFrame): A pandas dataframe

    Returns:
        rows (int): number of rows in input dataframe
        columns (int): number of columns in input dataframe
        row_lab (list): labels of rows
        col_lab (list): labels of columns
    """
    size = data.shape
    rows = size[0]
    columns = size[1]
    row_lab = data.index.tolist()
    col_lab = __make_cepnames(data.columns)

    with open('cep.dat', "w") as outfile:
        title = 'dummy title\n'
        outfile.write(title)

        data_format = '(I3,5(I3,F3.0))\n'
        outfile.write(data_format)

        outfile.write("%3d \n" % 5)

        for i in range(rows):
            number_written = 0
            row_number_written = False
            for j in range(columns):
                if data.iloc[i, j] != 0:
                    if not row_number_written:
                        row_number_written = True
                        outfile.write("%3d" % (i + 1))
                    outfile.write("%3d%3.0f" % ((j + 1), data.iloc[i, j]))
                    number_written += 1
                    if number_written == 5:
                        outfile.write('\n')
                        # outfile.write("%3d" % (i + 1))
                        number_written = 0
                        row_number_written = False
            if number_written != 0:
                outfile.write('\n')

        outfile.write("%3d \n" % 0)

        for j, m in enumerate(col_lab, 1):
            outfile.write(m + ['', '\n'][j % 10 == 0])
        outfile.write('\n')
        for i, m in enumerate(row_lab, 1):
            outfile.write(str(m).ljust(8) + ['', '\n'][i % 10 == 0])

    return rows, columns, row_lab, col_lab


def __make_cepnames(names, seconditem=False):
    """Abbreviates a botanical or zoological Latin name into an eight-character name.

    Cornell Ecology Programs (CEP) use eight-letter abbreviations for species
    and site names. In species, the names are formed by taking four first
    letters of the generic name and four first letters of the specific or
    subspecific epithet. In this function, the CEP name is made by taking the
    four first letters of the first element, and four first letters of the last
    (default) or the second element (with seconditem=True). If there was only
    one name element, it is abbreviated to eight letters. The returned names
    are made unique by adding numbers to duplicate names. These names may be
    practical to avoid congestion in ordination plots.

    This function has been fully taken from https://github.com/maurobio/cornpy.

    Args:
        names (list): a list of names to be formatted into CEP names
        seconditem (bool): Take always the second item of the original name to the
        abbreviated name instead of the last original item (default is False).

    Returns:
        cepnames (list): a list of CEP names
    """
    count = 1
    cepnames = list()
    for name in names:
        n = len(name.split()) if not seconditem else 2
        str1 = name.split()[0]
        try:
            str2 = name.split()[n - 1]
        except IndexError:
            str2 = str1
        cepstr = str1[0:4] + str2[0:4].rstrip('.') if str1 != str2 else str1[0:8]
        if cepstr not in cepnames:
            cepnames.append(cepstr)
        else:
            cepnames.append(cepstr[0:7] + str(count))
            count += 1
    return cepnames


def main():
    """Example code for detrended correspondence analysis.

    This functions loads a sample dataset, perform the DCA and plot the result as a biplot.
    """
    df = pd.read_csv('./data/gauch.csv', index_col=0)
    site_scores, species_scores, site_labels, species_labels = decorana(df)
    biplot(site_scores, species_scores, site_labels, species_labels)


if __name__ == '__main__':
    main()