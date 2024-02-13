# DPXRust

![Image of protein (3ZYZ) that has been run through DPXRust](example.png)
*Image of 3ZYZ after being run through DPXRust colorized by DPX value*

This is a rewrite of the DPX algorithm in Rust. DPX is an algorithm for computing a score that indicates to what degree an atom is buried in a protein structure by finding the distance between every atom and the nearest atom with an ASA/SASA above a threshold (Defaults to 150.0 Angstroms squared) [1].

## Usage

Output CSV Example:
```shell
./DPXRust --input-path my_protein.pdb --csv-output-path output.csv
```

Output PDB Example (When outputting as PDB DPX values are saved in B-factor column):
```shell
./DPXRust --input-path my_protein.pdb --pdb-output-path output.pdb
```

Since version 1.1.0 DPXRust automaticlly generates SASA values for input proteins. But if you want to provide your own custom SASA values then you can add the SASA value for each atom to the B-factor field of the input file.

Note: DPXRust also supports mmCIF files

Note: Make sure to split protein complexes into different files before running DPXRust.

For a full list of CLI arguemnts run DPXRust with `--help`

## Citations:

[1]: Pintar A., Carugo O., Pongor S. DPX: for the analysis of the protein core. Bioinformatics. 2003;19:313–314
