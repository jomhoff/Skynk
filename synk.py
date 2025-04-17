#!/usr/bin/env python3

"""
Skynk Syntenty Plotting Tool

This script combines BUSCO gene output from Compleasm and karyotype data from two species to generate
a synteny ideogram using RIdeogram.

USAGE EXAMPLE:
--------------
python skynk.py \
  --karyotype1 pfas_karyotype.txt \
  --karyotype2 tiliqua_karyotype.txt \
  --busco1 pfas_busco_format.tsv \
  --busco2 tiliqua_busco_format.tsv \
  --rep1 replacement_pfas.txt \
  --rep2 replacement_tiliqua.txt \
  --outdir output_directory \
  --cmap viridis \
  --plot

REQUIRED INPUTS:
----------------
--karyotype1       First species karyotype file (must have columns: Chr, Start, End, Species)
--karyotype2       Second species karyotype file (same format as above)
--busco1           BUSCO full table from species 1
--busco2           BUSCO full table from species 2
--rep1             Replacement map for species 1 chromosome names (ie. chr1	1)
--rep2             Replacement map for species 2 chromosome names
--outdir           Directory for output files

OPTIONAL:
---------
--cmap             Colormap for chromosome gradient (default: plasma)
--plot             Include this flag to generate an ideogram PNG/SVG with RIdeogram
--rscript_path     Path to output R script (default: plot_ideogram.R)
--karyo_size       Label text size (default: 12)
--karyo_color      Label text color (default: black)

Output Files:
-------------
- chr_color_map.txt      : Map of chromosome names to colors
- color_replace.txt      : Map of numeric chromosome order to colors
- merged_busco.txt       : Merged BUSCO hits shared between both species
- final_synteny.txt      : Synteny links colored by species 1 chromosome
- dual_karyotype.txt     : Combined karyotype file with label and fill info
- chromosome.png/.svg    : RIdeogram synteny plot (if --plot is used)

Author: Jon J. Hoffman

Copyright 2025 Jon J. Hoffman

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


# Python dependencies
# pip install pandas matplotlib

# Or with conda
# conda create -n synteny_env python=3.9 pandas matplotlib
# conda activate synteny_env

import argparse
import os
import sys
import csv
import pandas as pd
from collections import OrderedDict
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import subprocess

# === VISUALS === #
ERROR_ART = """
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⣿⣷⠀⠀⠀⠀⣸⣶⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣼⣿⣿⡞⣿⣷⣮⣻⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⣿⣿⣿⣿⣾⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⡝⢿⣿⣿⣿⣿⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣤⡀⠀⠀⠀⠀⠀⠀⠻⣿⣿⣿⠸⣸⣻⣏⣿⣿⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⣿⣿⡿⡀⠀⠀⠀⠀⠀⣾⡞⡝⣿⢿⣿⣿⣿⣿⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠩⣾⣿⣶⢦⣤⣀⠸⠻⢭⣥⡻⣧⠀⡙⠛⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⣄⢠⣴⣾⣿⣿⣿⣏⣶⣾⡽⣿⣷⣟⣿⣿⣿⣻⣷⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⣀⣀⣀⠀⠀⠀⠸⣿⡿⠘⠻⢿⣿⣿⠟⠛⠿⠿⠃⢍⣿⣿⢸⣿⣿⣿⡽⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⣰⣟⠛⠛⢿⣿⣦⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣿⣜⢿⣿⡿⡷⡿⣼⣶⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⢰⣿⠃⠀⠀⠀⠈⢿⣿⣧⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣷⣯⣾⣿⡀⠀⠙⠻⢿⣶⣄⠀⠀⠀⠀⠀⠀⠀
⢸⣿⠀⠀⠀⠀⠀⠀⢻⣿⣷⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣧⡀⠀⠀⠀⠙⢿⣧⡀⠀⠀⠀⠀⠀
⢸⣿⡀⠀⠀⠀⠀⠀⠀⢻⣿⣿⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⣬⣽⣿⣿⢟⣛⣳⠀⠀⠀⠀⠀⠹⣿⣆⠀⠀⠀⠀
⠀⣿⣇⠀⠀⠀⠀⠀⠀⠈⣿⣿⣿⣷⡀⠀⠀⠀⠀⠀⠀⠀⣴⣿⣿⣿⣿⣷⢻⣾⣿⣿⣷⡽⣄⠀⠀⢀⣾⣿⣷⣄⠀⠀
⠀⠘⣿⣆⠀⠀⠀⠀⠀⠀⠘⣿⣿⣿⣿⣷⣄⡀⠀⠀⢀⣾⣿⣿⣿⣿⣿⣿⡇⣿⣿⣿⣿⣿⢹⣦⠀⢸⣇⠀⠹⣏⢧⡀
⠀⠀⠹⣿⣷⡀⠀⠀⠀⠀⠀⠘⣿⣿⣿⣿⣿⣿⣿⡆⣿⣿⣿⣿⣿⣿⣿⣿⣧⣿⣿⣿⣿⣿⢸⣿⡄⠈⠛⠀⣶⠟⠼⠇
⠀⠀⠀⠹⣿⣿⣷⣤⡀⠀⠀⠀⠘⢿⣿⣿⣿⣿⣿⢸⣿⣿⣿⣿⣿⣿⣿⡿⣼⣿⣿⣿⣿⡿⣾⣿⠁⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠙⣿⣿⣿⣿⣶⣄⠀⠀⠈⠻⣿⣿⣿⣿⢸⣿⣿⣿⣿⣿⣿⡿⣱⣿⣿⣿⣿⢟⣼⣿⠏⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠈⢻⣿⣿⣿⣿⣧⡀⠀⠀⠈⠻⢿⣿⢸⣿⣿⣿⡿⢟⣫⣾⣿⣿⠿⣛⣵⣿⡿⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠙⠿⣿⣿⣿⡇⠀⠀⠀⠀⠀⢈⣾⣿⡟⠙⠚⠛⠛⠋⠉⠀⠘⣿⣿⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠛⠁⠀⠀⠀⠀⢀⣾⣿⡟⠀⠀⠀⠀⠀⠀⠀⠀⢰⣿⡿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣾⣿⠏⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⣿⣷⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣿⡿⡏⠀⠀⠀⠀⠀⠀⠀⠀⢠⣾⣯⢻⣿⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⣿⣿⣧⡀⠀⠀⠀⠀⠀⠀⠀⠈⠛⠋⠘⠻⣿⣿⣷⣶⣒⣒⢢⡄⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⣿⣿⣿⡿⣏⣃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠻⠿⠿⠟⠈⠁⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⣿⡿⠿⠿⠿⣿⣿⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
"""

SUCCESS_ART = """
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⡀⠀⠄⣀⡀⡀⠤⠐⣢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⠠⢤⠔⠈⠀⠀⠀⠀⠀⠀⠀⠁⠀⣾⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢳⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢛⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢿⠏⠀⣀⣀⡀⠀⠀⠀⠀⢀⠀⡔⢻⣦⠀⢃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⣀⡀⠀⠀⠀⠀⣀⡀⠀⠀⠀⠀⠀⠀⠸⠀⠰⠛⣇⣹⡜⡄⠀⠀⠸⢠⣿⣿⣀⠇⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⢀⠔⡉⠤⠐⠒⠒⠒⠂⠠⠬⣁⠒⠠⢄⡀⠀⢠⠀⠐⠠⠿⢿⡇⠀⠀⠀⠀⠈⠛⠉⠀⡀⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⡐⡡⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠐⠂⠄⡁⠒⠱⢤⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠤⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠤⠄⠒⣒⣀⣴
⠰⠰⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠐⠢⠌⣉⣶⣶⣦⣄⣀⣠⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠠⠂⢁⣠⣴⣾⣿⣿⡟⠁
⡇⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠀⠉⠛⠿⠿⠛⠆⠤⠄⠀⣀⣀⣀⣀⣀⡀⠔⢁⣤⣾⣿⣿⣿⡿⠟⠁⠀⠀
⡇⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠁⢀⠀⠀⠄⠀⡀⠀⠁⠐⠒⠂⠠⠤⢤⣤⣤⣶⣾⣿⠿⠿⠛⠋⠁⠀⠀⠀⠀⠀
⢇⢃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠁⠀⠀⠈⡄⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠘⡈⢆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠁⠀⠀⠉⢂⠀⢱⡀⡰⢠⡃⣸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠐⢄⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣿⠀⠀⠀⠀⢸⠀⠀⠈⠀⠀⠉⡉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠑⢄⡁⠂⠤⢀⣀⠀⠀⠀⠀⣀⣀⣼⣿⠀⠀⠀⠀⣮⣤⣶⣶⣦⠋⢀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠈⠑⠂⠤⠄⠀⠀⠠⠾⠿⠿⢻⣿⠀⠀⢀⣴⣿⣿⣿⡿⠁⡀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣶⡿⢃⠤⠐⠁⠀⢰⣾⠟⡀⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⡏⠀⠆⠀⠀⠀⠀⢸⡟⠀⢡⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣿⠀⢰⠀⠀⠀⠀⠀⠀⣇⠀⠈⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡟⠀⠸⠀⠀⠀⠀⠀⠀⢹⠀⠀⢁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⡇⠀⡀⠀⠀⠀⠀⠀⠀⢸⣇⠀⠘⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⣿⡇⠀⡇⠀⠀⠀⠀⠀⠀⢸⣿⡀⠀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⡿⠁⣆⠁⠀⠀⠀⠀⠀⠀⠀⣿⣷⢠⡸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠉⠁⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
"""

def fail(message):
    print(ERROR_ART)
    print(f"\n\U0001F4A5 {message}")
    sys.exit(1)

def succeed(message):
    print(SUCCESS_ART)
    print(f"\n\u2705 {message}")

def log(message, icon="\U0001F527"):
    print(f"{icon} {message}")

def ensure_outdir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        log(f"Created output directory: {path}", icon="\U0001F4C1")

# === COLOR MAP ===
ROMAN_NUMERAL_MAP = {
    'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8,
    'IX': 9, 'X': 10, 'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15,
    'XVI': 16, 'XVII': 17, 'XVIII': 18, 'XIX': 19, 'XX': 20
}

def try_parse_chr(chr_label):
    try:
        return int(chr_label)
    except ValueError:
        return ROMAN_NUMERAL_MAP.get(chr_label.upper(), chr_label)

#generate color gradient
def generate_color_gradient(num_colors, cmap_name='plasma'):
    cmap = cm.get_cmap(cmap_name, num_colors)
    return [mcolors.to_hex(cmap(i)).replace('#', '') for i in range(num_colors)]

#Generate Color map
def create_chr_color_map(karyo_file, output_file, replace_output_file, cmap_name='plasma'):
    try:
        df = pd.read_csv(karyo_file, sep="\t")
        if 'Chr' not in df.columns:
            fail("Karyotype file must contain a 'Chr' column.")
        chromosomes = df['Chr'].astype(str).unique()
        chromosomes_sorted = sorted(chromosomes, key=try_parse_chr)
        colors = generate_color_gradient(len(chromosomes_sorted), cmap_name)

        with open(output_file, 'w') as out, open(replace_output_file, 'w') as rep:
            for i, (chr_label, color) in enumerate(zip(chromosomes_sorted, colors)):
                out.write(f"{chr_label}\t{color}\n")
                rep.write(f"{i+1}\t{color}\n")
        log(f"Wrote chromosome-color map to: {output_file}", icon="\U0001F3A8")
    except Exception as e:
        fail(f"Failed to generate color map: {e}")

def filter_non_integer_chrs(filepath):
    try:
        cleaned = []
        with open(filepath, "r") as f:
            for i, line in enumerate(f):
                parts = line.strip().split("\t")
                if i == 0 or (parts[0].isdigit() and parts[3].isdigit()):
                    cleaned.append(line)
                else:
                    log(f"Removing line {i+1} with non-integer Chr: {line.strip()}", icon="\U0001F9F9")
        with open(filepath, "w") as f:
            f.writelines(cleaned)
        log("Removed rows with non-integer chromosome labels.", icon="\u2705")
    except Exception as e:
        fail(f"Failed to filter non-integer chromosome labels: {e}")

#combine karyotypes
def augment_karyotype(karyo_file, color_map, fill_all=None, size=None, color=None):
    try:
        df = pd.read_csv(karyo_file, sep="\t")

        if fill_all:
            df['fill'] = fill_all
        else:
            chr_color_map = pd.read_csv(color_map, sep="\t", header=None, names=["Chr", "fill"])
            df = df.merge(chr_color_map, on="Chr", how="left")

            # Safely resolve fill columns
            fill_x = df['fill_x'] if 'fill_x' in df.columns else pd.Series([None] * len(df))
            fill_y = df['fill_y'] if 'fill_y' in df.columns else pd.Series([None] * len(df))
            fill   = df['fill']    if 'fill' in df.columns    else pd.Series([None] * len(df))

            df['fill'] = fill_x.combine_first(fill_y).combine_first(fill)
            df.drop(columns=[col for col in ['fill_x', 'fill_y'] if col in df.columns], inplace=True)

        df['size'] = size
        df['color'] = color

        # Enforce final column order expected by RIdeogram
        expected_order = ['Chr', 'Start', 'End', 'fill', 'species', 'size', 'color']
        df = df[[col for col in expected_order if col in df.columns]]

        df.to_csv(karyo_file, sep="\t", index=False)
        log(f"🧰 Augmented karyotype file: {karyo_file}", icon="🧰")

    except Exception as e:
        fail(f"Error augmenting {karyo_file}: {e}")

# === BUSCO MERGE ===
def merge_busco(busco1_path, busco2_path, output_path):
    try:
        def read_busco(path, label):
            data = {}
            with open(path) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    if row.get("Status", "").strip() != "Complete":
                        continue
                    bid = row["# Busco id"].strip().lower()
                    data.setdefault(bid, {})[label] = row
            return data

        b1 = read_busco(busco1_path, "1")
        b2 = read_busco(busco2_path, "2")
        shared = [bid for bid in b1 if bid in b2]
        print(f"🧬 Shared complete BUSCOs: {len(shared)}")

        with open(output_path, "w") as out:
            out.write("Species_1\tStart_1\tEnd_1\tSpecies_2\tStart_2\tEnd_2\tfill\n")
            for bid in shared:
                try:
                    row1 = b1[bid]["1"]
                    row2 = b2[bid]["2"]
                    chr1 = row1.get("Contig") or row1.get("Sequence")
                    chr2 = row2.get("Contig") or row2.get("Sequence")
                    start1 = int(float(row1["Gene Start"]))
                    end1 = int(float(row1["Gene End"]))
                    start2 = int(float(row2["Gene Start"]))
                    end2 = int(float(row2["Gene End"]))
                    if not chr1 or not chr2:
                        continue
                    out.write(f"{chr1}\t{start1}\t{end1}\t{chr2}\t{start2}\t{end2}\tplaceholder\n")
                except:
                    continue
    except Exception as e:
        fail(f"Failed to merge BUSCO files: {e}")


def apply_chr_replacements(synteny_file, rep1_path, rep2_path):
    try:
        def load_map(path):
            m = {}
            with open(path) as f:
                for line in f:
                    old, new = line.strip().split("\t")
                    m[old] = new
            return m

        rep1 = load_map(rep1_path)
        rep2 = load_map(rep2_path)

        out_lines = []
        with open(synteny_file) as f:
            header = f.readline().strip()
            out_lines.append(header)
            for i, line in enumerate(f, 2):
                parts = line.strip().split("\t")
                if len(parts) != 7:
                    continue
                parts[0] = rep1.get(parts[0], parts[0])
                parts[3] = rep2.get(parts[3], parts[3])
                out_lines.append("\t".join(parts))

        with open(synteny_file, "w") as f:
            for line in out_lines:
                f.write(line + "\n")
        log("Chromosome names replaced using provided maps.", icon="\U0001F504")
    except Exception as e:
        fail(f"Failed to apply chromosome name replacements: {e}")

def replace_fill(input_path, color_map_path, output_path):
    try:
        color_dict = {}
        with open(color_map_path) as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        color_dict[parts[0]] = parts[1]

        with open(input_path) as infile, open(output_path, "w") as outfile:
            for i, line in enumerate(infile):
                if not line.strip():
                    continue
                parts = line.strip().split("\t")
                if i == 0:
                    outfile.write(line)
                    continue
                if len(parts) != 7:
                    continue
                key = parts[0]
                fill = color_dict.get(key)
                if fill:
                    parts[-1] = fill
                    outfile.write("\t".join(parts) + "\n")
    except Exception as e:
        fail(f"Failed to replace fill values: {e}")

def write_and_run_rscript(karyotype_path, synteny_path, working_dir, script_path):
    r_script = f"""
if (!requireNamespace("RIdeogram", quietly=TRUE)) {{
  if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
  devtools::install_github("TickingClock1992/RIdeogram")
}}
library(RIdeogram)
if (!requireNamespace("rsvg", quietly=TRUE)) install.packages("rsvg")
library(rsvg)

setwd("{working_dir}")
print(getwd())
print(list.files())

kary <- read.table("{os.path.basename(karyotype_path)}", header=TRUE)
synt <- read.table("final_synteny.txt", header=TRUE, colClasses = c("numeric", "integer", "integer", "numeric", "integer", "integer", "character"))

ideogram(karyotype=kary, synteny=synt)
rsvg_png("chromosome.svg", "chromosome.png", width=1000)
"""
    with open(script_path, "w") as f:
        f.write(r_script)

    result = subprocess.run(["Rscript", script_path], capture_output=True)
    if result.returncode != 0 or not os.path.exists(os.path.join(working_dir, "chromosome.png")):
        fail(f"R plotting failed: {result.stderr.decode()}")
    else:
        log("RIdeogram plot created successfully!", icon="\U0001F3A8")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--karyotype1", required=True)
    parser.add_argument("--karyotype2", required=True)
    parser.add_argument("--busco1", required=True)
    parser.add_argument("--busco2", required=True)
    parser.add_argument("--rep1", required=True)
    parser.add_argument("--rep2", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--cmap", default="plasma")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--rscript_path", default="plot_ideogram.R")
    parser.add_argument("--karyo_size", default="12")
    parser.add_argument("--karyo_color", default="black")
    args = parser.parse_args()

    ensure_outdir(args.outdir)

    chr_color_map = os.path.join(args.outdir, "chr_color_map.txt")
    color_replace_map = os.path.join(args.outdir, "color_replace.txt")
    combined_karyotype = os.path.join(args.outdir, "combined_karyotype.txt")
    merged_busco = os.path.join(args.outdir, "merged_busco.txt")
    final_synteny = os.path.join(args.outdir, "final_synteny.txt")

    create_chr_color_map(args.karyotype1, chr_color_map, color_replace_map, cmap_name=args.cmap)
    augment_karyotype(args.karyotype1, chr_color_map, size=args.karyo_size, color=args.karyo_color)
    augment_karyotype(args.karyotype2, chr_color_map, fill_all="cccccc", size=args.karyo_size, color=args.karyo_color)


    df1 = pd.read_csv(args.karyotype1, sep="\t")
    df2 = pd.read_csv(args.karyotype2, sep="\t")
    combined = pd.concat([df1, df2], ignore_index=True)
    combined_karyo_path = os.path.join(args.outdir, "dual_karyotype.txt")
    combined.to_csv(combined_karyo_path, sep="\t", index=False)
    log(f"Wrote combined karyotype: {combined_karyo_path}", icon="\U0001F91D")

    merge_busco(args.busco1, args.busco2, merged_busco)
    apply_chr_replacements(merged_busco, args.rep1, args.rep2)
    filter_non_integer_chrs(merged_busco)
    replace_fill(merged_busco, color_replace_map, final_synteny)

    if args.plot:
        write_and_run_rscript(
            karyotype_path=combined_karyo_path,
            synteny_path=final_synteny,
            working_dir=args.outdir,
            script_path=os.path.join(args.outdir, args.rscript_path)
        )

    succeed(f"Pipeline complete! Final synteny file: {final_synteny}")

if __name__ == "__main__":
    main()
