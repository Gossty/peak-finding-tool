# Peak Finding

This repository is used to find peaks in a provided tag directory. It takes as input a tag directory containing tag files, along with a control directory for peak finding. The code performs various filtering and analysis steps to identify peaks and outputs the results in a BED file format.


## Usage

```
python peakFinding.py [tag_directory] [control] [options]
```

Optional parameters:
```
    [-style <factor/histone>] [-o <output_directory>] [-poisson <threshold>] [-fold <fold_change>] [-fragLen <window_length>]
```

- `tag_directory`: Path to the tag directory containing the tag files for analysis.
- `control`: Path to the control directory for peak finding.
- `-style <factor/histone>` (optional): Specify the style. Use `<factor>` for TFs or `<histone>` for histone modifications.
- `-o <output_directory>` (optional): Specify the output directory for the results. If not provided, the results will be saved to the current directory.
- `-poisson <threshold>` (optional): Manually set the threshold for the Poisson filter. Default is 1e-60.
- `-fold <fold_change>` (optional): Manually set the fold change for peak detection. Default is 4.
- `-fragLen <window_length>` (optional): Manually specify the length of the window. Default is 75.

## Workflow

1. The code reads the command-line arguments to gather the necessary input for analysis.
2. It processes the tag directory and control directory, filtering for necessary columns and storing the data in dictionaries.
3. The code calculates the tag counts for each window for both the sample and control data, saving these to a dictionary where the key is the position and the value is the count associated with that position.
4. Filtering steps are applied to identify potential peaks:
   - Fold Change Filtering: The code filters based on the fold change, where a fold value of 4 or greater indicates a binding site.
   - Maximum Count Filtering: The code selects the position with the maximum tag count within each window of width WINDOW_LENGTH * 2, so as to not double count potential peaks.
   - Poisson Filtering: The code applies a Poisson filter using the expected value and a threshold to further refine the peak selection. It calculates the probability of seeing the observed number of tags at a given window using expected tags per window obtained in the control files. Expected value λ= (WL * #iput_tags) / GL, where WL is the fragment length estimate from makeTagDirectory, #input_tags is the number of tags in the input genome, and GL is the estimated length of the genome (2e9)
5. The final filtered peaks are saved in a BED file format named `example.bed` in the specified output directory (or the current directory if not provided).

Note: The code includes print statements for displaying various statistics and intermediate results during the analysis process.

## Example

```
python peakFinding.py ./tag_directory ./control
```

This example command runs the peak finding analysis on the tag files in the `tag_directory` using the `control` directory. 

## Data

`example_files`: directory containing tag directories for the first 16 chromosomes for input and Klf4 transcription factor. This data should return an empty `bed` file because it is noise data.
`tests`: directory containing the test data of first 17 chromosomes for input and Klf4 transcription factor. This data should return a `bed	` file indicating all the peaks.

The dataset is originally from this paper:
Chronis C, Fiziev P, Papp B, Butz S, Bonora G, Sabri S, Ernst J, Plath K. Cooperative Binding of Transcription Factors Orchestrates Reprogramming. Cell. 2017 Jan 26;168(3):442-459.e20. doi: 10.1016/j.cell.2016.12.016. Epub 2017 Jan 19. PMID: 28111071; PMCID: PMC5302508.
https://www.ncbi.nlm.nih.gov/pubmed/28111071

##Credit:
Our current implementation heavily relies on documentation of the following libraries:

MACS2
HOMER
Thanks to Prof. Gymrek and TA Ryan at UCSD for guiding us through the project.

## Contact

For any issues or questions, please contact Stephen Golzari, Yashwin Madakamutil, or Yang Han at sgolzari@ucsd.edu, yamadamutil@ucsd.edu, or yah015@ucsd.edu .
