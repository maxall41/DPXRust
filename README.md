# DPXRust

![Image of protein that has been run through DPXRust](example.png)
*Image of 3ZYZ after being run through DPXRust colorized by DPX value*

This is a rewrite of the DPX algorithm in Rust. DPX is an algorithm for computing a score that indicates to what degree an atom is buried in a protein structure by finding the distance between every atom and the nearest atom with an ASA/SASA above a threshold (Defaults to 10.0 Angstroms squared) [1].

## Usage

Make sure that the SASA/ASA value for each atom is in the B-factor field before running DPXRust.

Output CSV Example:
```shell
./DPXRust --input-path my_protein.pdb --csv-output-path output.csv
```

Output PDB Example (When outputting as PDB DPX values are saved in B-factor column):
```shell
./DPXRust --input-path my_protein.pdb --pdb-output-path output.csv
```

Note: DPXRust also supports mmCIF files

## Citations:

[1]: Pintar A., Carugo O., Pongor S. DPX: for the analysis of the protein core. Bioinformatics. 2003;19:313–314
