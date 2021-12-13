"""convert values from HELIOS opacity data file from CGS to SI"""

from alfrodull_input_tools import convert_helios_opacity_to_SI_units


import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='Alfrodull unit converter')

parser.add_argument("source_file", action="store", help = "input file to read data from, in CGS units")
parser.add_argument("output_file", action="store", help = "output file to write data in SI units")

args = parser.parse_args()

source_file = Path(args.source_file)
output_file = Path(args.output_file)

if not source_file.exists():
    print(f"Error: {source_file} input does not exist")
    exit(-1)

if output_file.exists():
    print(f"Error: {output_file} output already exist")
    exit(-1)

convert_helios_opacity_to_SI_units(source_file,output_file)

print("Conversion done")
