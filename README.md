# DPXRust

This is a rewrite of the DPX algorithm in Rust. DPX is an algorthim for computing the depth of residues in a protein structure by finding the distance between the current residue and the nearest residue with a SASA above 0.0 [1].

## Usage

Example:
```shell
./DPXRust --input-path my_protein.pdb --output-path output.csv
```
## Citations:

[1]: Pintar A., Carugo O., Pongor S. DPX: for the analysis of the protein core. Bioinformatics. 2003;19:313–314