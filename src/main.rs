use std::cmp::Ordering;
use std::fs::File;
use pdbtbx::*;
use polars::prelude::*;

#[derive(Debug)]
struct ASAData {
    asa: f64,
    position: (f64,f64,f64),
    residue_id: i32
}

#[derive(Debug)]
struct DPXData {
    dpx: f64,
    residue_id: i32
}

fn distance(p1: &(f64, f64, f64), p2: &(f64, f64, f64)) -> f64 {
    let (x1, y1, z1) = p1;
    let (x2, y2, z2) = p2;

    let dx = x1 - x2;
    let dy = y1 - y2;
    let dz = z1 - z2;

    (dx*dx + dy*dy + dz*dz).sqrt()
}

fn sort_by_distance(data: &mut Vec<ASAData>, reference_position: (f64, f64, f64)) {
    data.sort_unstable_by(|a, b| {
        let da = distance(&a.position, &reference_position);
        let db = distance(&b.position, &reference_position);

        if da < db {
            Ordering::Less
        } else if da > db {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    });
}

fn calculate_atomic_info_for_residue(residue: &Residue) -> (f64,(f64,f64,f64)) {
    let mut average_asa: f64 = 0.0;
    let mut average_residue_position = (0.0,0.0,0.0);
    for atom in residue.atoms() {
        average_asa += atom.b_factor();
        average_residue_position.0 = average_residue_position.0 + atom.x();
        average_residue_position.1 = average_residue_position.1 + atom.y();
        average_residue_position.2 = average_residue_position.2 + atom.z();
    }
    average_asa = average_asa / residue.atom_count() as f64;
    average_residue_position.0 = average_residue_position.0 / residue.atom_count() as f64;
    average_residue_position.1 = average_residue_position.1 / residue.atom_count() as f64;
    average_residue_position.2 = average_residue_position.2 / residue.atom_count() as f64;
    return (format!("{:.2}", average_asa).parse().unwrap(), average_residue_position);
}
fn main() {
    const EXPOSED_ASA_THRESHOLD: f64 = 0.0; // Units is angstroms squared
    println!("WARNING: This program expects the B-factor field of the input PDB file to be the ASA/SASA of each residue. IF THIS IS NOT THE CASE THIS PROGRAM WILL NOT WORK.");
    let (mut pdb, _errors) = pdbtbx::open(
        "./example.pdb",
        StrictnessLevel::Medium
    ).unwrap();

    pdb.remove_atoms_by(|atom| atom.element() == Some(&Element::H)); // Remove all H atoms
    use std::time::Instant;
    let now = Instant::now();
    // First loop to build list of ASA values for each residue
    let mut asa_values : Vec<ASAData> = Vec::new();
    for residue in pdb.residues() { // Iterate over all residues in the structure
        let (average_asa, position) = calculate_atomic_info_for_residue(residue);
        let data = ASAData {
            asa: average_asa,
            residue_id: residue.id().0 as i32,
            position
        };
        asa_values.push(data);
    }
    // Second loop to generate DPX values for each residue
    let mut dpx_residue_ids: Vec<i32> = vec![];
    let mut dpx_values: Vec<f64> = vec![];
    for residue in pdb.residues() {
        let (average_asa,position) = calculate_atomic_info_for_residue(residue);
        if average_asa < EXPOSED_ASA_THRESHOLD {
            // Atom is not exposed
            sort_by_distance(&mut asa_values,position);
            for asa_value in &asa_values {
                if asa_value.residue_id == residue.id().0 as i32 {
                    continue
                }
                if asa_value.asa > EXPOSED_ASA_THRESHOLD {
                    let dpx = distance(&position, &asa_value.position);
                    dpx_residue_ids.push(residue.id().0 as i32);
                    dpx_values.push(dpx);
                    break
                }
            }
        } else {
            // Atom is exposed
            dpx_residue_ids.push(residue.id().0 as i32);
            dpx_values.push(0.0);
        }
    }
    let elapsed = now.elapsed();
    let mut df: DataFrame = df!(
        "Residue ID" => &dpx_residue_ids,
        "DPX" => &dpx_values
    ).unwrap();
    println!("{}", df);
    let mut file = File::create("./output.csv").expect("could not create file");
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b',')
        .finish(&mut df).unwrap();
}
