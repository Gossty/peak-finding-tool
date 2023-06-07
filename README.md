# Peak Finding
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

This repository is used to find peaks in a provided tag directory. It takes as input a tag directory containing tag files, and it optionally takes a control directory for peak finding as well as other parameters. The code performs various filtering and analysis steps to identify peaks and outputs the results in a BED file format and a txt file with statistical information.

## Required Packages
* numpy
* pandas
* scipy
* argparse

## Install Instructions
Once the required libraries are installed, you can install `peakFinding` by running the following command inside of `peak-finding-tool` directory:

`python setup.py install --user`

If install was successful, typing `peakFinding --help` will print the usage information of the tool. After that you would be able to use the tool.

Alternatively, you can install this tool by running the following commands inside of `peak-finding-tool` directory:

`pip install wheel sdist `

`python setup.py bdist_wheel sdist`

`pip install .`

If the following commands give you a warning that the installation was not on PATH, you would need to add the directory to which it installed into your shell's path (on mac for zsh it would be `path+=<the_path_where_saved>`; on bash `export PATH=$PATH:<the_path_where_saved`, e.g. `export PATH=$PATH:$HOME/.local/bin`)



**Note:**

If you are unable to install the `peakFinding` tool, you can take another way to run our code: (please make sure that the paths are correct on your device)

Here is an example of running it in the terminal
```
git clone https://github.com/Gossty/peak-finding-tool.git

cd peakFinding 

python peakFinding.py [tag_directory] [options]
```

**If you run the code using this way (through .py file), make sure to use the output directory inside of `./peakFinding/` when running with `-o` option (i.e. create a new directory inside of `./peakFinding/`).**

## Usage

First of all you need to obtain a tag_directory before using the tool. You can do that by using a command line from Homer:

```
makeTagDirectory <output directory> <input BAM file> [options]
```

You can use sample test data inside of `./tests/` and `./example-files/` directory of GitHub repository.

```
peakFinding [tag_directory] [options]
```

Optional parameters:
```
    [-control <control_directory>] [-o <output_directory>] [-poisson <threshold>] [-fold <fold_change>] [-fragLen <window_length>]
```

- `tag_directory`: Path to the tag directory containing the tag files for analysis.
- `-control`: (optional) Path to the control directory for peak finding.
- `-o <output_directory>` (optional): Specify the output directory for the results. If not provided, the results will be saved to the current directory.
- `-poisson <threshold>` (optional): Manually set the threshold for the Poisson filter. Default is 1e-4.
- `-fold <fold_change>` (optional): Manually set the fold change for peak detection. Default is 4.
- `-fragLen <window_length>` (optional): Manually specify the length of the window. Default is 75.
- `-L <local_window>` (optional): Manually specify the length of the window for local filter. Default is 10000.


## Workflow

1. The code reads the command-line arguments to gather the necessary input for analysis.
2. It processes the tag directory, filtering for necessary columns and storing the tag data in dataframes.
3. The code calculates the tag counts for each window for the sample data, saving these to a dictionary where the key is the position and the value is the count associated with that position.
4. Next are filtering steps described below.

### Filtering steps follow two separate workflows
#### **Without control/input directory**:


The following filtering steps are applied:
   - Maximum Count Filtering: The code selects the position with the maximum tag count within each window of width `WINDOW_LENGTH * 2`, so as to not double count potential peaks.
   - Local Density Filtering: The code filters based on the fold change, of the current density to the local density of window with width `LOCAL_WINDOW`; filtering occurs based on the fold value specified with `FOLD_VALUE`.
   - Poisson Filtering: The code applies a Poisson filter using the expected value and a threshold to further specify the peak selection. It calculates the probability of seeing the observed number of tags at a given window using expected tags per window obtained in the local files and filters the peaks through the specified `THRESHOLD`. Expected value `exp = (WINDOW_LENGTH * #sample_tags) / GENOME_LENGTH`. 
#### **With input directory**: 

Filtering steps are applied to identify potential peaks:
   - Maximum Count Filtering: The code selects the position with the maximum tag count within each window of width WINDOW_LENGTH * 2, so as to not double count potential peaks.
   - Fold Change Filtering: The code filters based on the fold change, where a fold value of `FOLD_VALUE` or greater indicates a binding site.
   - Local Density Filtering: The code filters based on the fold change, of the current density to the local density of window with width `LOCAL_WINDOW`; filtering occurs based on the fold value specified with `FOLD_VALUE`.
   - Poisson Filtering: The code applies a Poisson filter using the expected value and a threshold to further specify the peak selection. It calculates the probability of seeing the observed number of tags at a given window using expected tags per window obtained in the local files and filters the peaks through the specified `THRESHOLD`. Expected value `exp = (WINDOW_LENGTH * #input_tags) / GENOME_LENGTH`. 
#
5. The final filtered peaks are saved in a BED file format named `peaks.bed` in the specified output directory (or the current directory if not provided). Additionaly outputs `peaks.txt` file with statistical data.

Note: The code includes print statements for displaying various intermediate results during the analysis process.

## Example

```
peakFinding ./tag_directory -control ./control_directory
```

This example command runs the peak finding analysis on the tag files in the `tag_directory` using the `control_directory`. 

## Data

`example-files`: directory containing tag directories for the first 16 chromosomes for input and Klf4 transcription factor. This data should return an empty `bed` file because it is noise data.

`tests`: directory containing the test data of all chromosomes for input and Klf4 transcription factors. This data should return a `peaks.bed` file indicating all the peaks and additional statistical information `peaks.txt`.

The dataset is originally from this paper:
Chronis C, Fiziev P, Papp B, Butz S, Bonora G, Sabri S, Ernst J, Plath K. Cooperative Binding of Transcription Factors Orchestrates Reprogramming. Cell. 2017 Jan 26;168(3):442-459.e20. doi: 10.1016/j.cell.2016.12.016. Epub 2017 Jan 19. PMID: 28111071; PMCID: PMC5302508.
https://www.ncbi.nlm.nih.gov/pubmed/28111071

## Credit
Our current implementation heavily relies on documentation of the following libraries:

- MACS2
- HOMER

Thanks to Prof. Gymrek and TA Ryan at UCSD for guiding us through the project.

## Contact

For any issues or questions, please contact Stephen Golzari, Yashwin Madakamutil, or Yang Han at sgolzari@ucsd.edu, yamadamutil@ucsd.edu, or yah015@ucsd.edu .
