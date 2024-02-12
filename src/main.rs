use clap::Parser;
use colored::Colorize;
use pdbtbx::*;
use polars::prelude::*;
use std::cmp::Ordering;
use std::fs::File;
use std::path::Path;
use std::process::exit;

#[derive(Debug)]
struct ASAData {
    asa: f64,
    position: (f64, f64, f64),
    residue_id: i32,
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the input PDB file
    #[arg(short, long)]
    input_path: String,

    /// Optional: Path for the CSV to be created
    #[arg(short, long)]
    csv_output_path: Option<String>,

    /// Optional: Path for the PDB to be created - Will not be created if path is not specified
    #[arg(short, long)]
    pdb_output_path: Option<String>,

    ///  Exposed ASA Threshold - Defaults to 10 angstroms squared (absolute)
    #[arg(short, long)]
    exposed_asa_threshold: Option<f64>,
}

// Calculate euclidean distance between points
fn distance(p1: &(f64, f64, f64), p2: &(f64, f64, f64)) -> f64 {
    let (x1, y1, z1) = p1;
    let (x2, y2, z2) = p2;

    let dx = x1 - x2;
    let dy = y1 - y2;
    let dz = z1 - z2;

    (dx * dx + dy * dy + dz * dz).sqrt()
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

fn calculate_atomic_info_for_residue(residue: &Residue) -> (f64, (f64, f64, f64)) {
    let mut average_asa: f64 = 0.0;
    let mut average_residue_position = (0.0, 0.0, 0.0);
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
    return (
        format!("{:.2}", average_asa).parse().unwrap(),
        average_residue_position,
    );
}
fn main() {
    let args = Args::parse();

    let mut exposed_asa_threshold = 10.0;

    if let Some(in_exposed_asa_threshold) = args.exposed_asa_threshold {
        exposed_asa_threshold = in_exposed_asa_threshold;
    }

    if Path::new(&args.input_path).exists() == false {
        println!(
            "{}",
            format!("ERROR: Input file ({}) does not exist", &args.input_path).red()
        );
        exit(1);
    }
    println!("{}", "WARNING: This program expects the B-factor field of the input PDB file to be the ASA/SASA of each residue. IF THIS IS NOT THE CASE THIS PROGRAM WILL NOT WORK.".yellow());
    let (mut pdb, _errors) = pdbtbx::open(args.input_path, StrictnessLevel::Loose).unwrap();
    if args.csv_output_path.is_none() && args.pdb_output_path.is_none() {
        println!(
            "{}",
            "You must specify --pdb-output-path and or --csv-output-path".red()
        );
        exit(1);
    }

    pdb.remove_atoms_by(|atom| atom.element() == Some(&Element::H)); // Remove all H atoms

    // First loop to build list of ASA values for each residue
    let mut asa_values: Vec<ASAData> = Vec::new();
    for residue in pdb.residues() {
        // Iterate over all residues in the structure
        let (average_asa, position) = calculate_atomic_info_for_residue(residue);
        let data = ASAData {
            asa: average_asa,
            residue_id: residue.id().0 as i32,
            position,
        };
        asa_values.push(data);
    }
    // Second loop to generate DPX values for each residue
    let mut dpx_residue_ids: Vec<i32> = vec![];
    let mut dpx_values: Vec<f64> = vec![];
    for residue in pdb.residues_mut() {
        let (average_asa, position) = calculate_atomic_info_for_residue(residue);
        if average_asa <= exposed_asa_threshold {
            // Atom is not exposed
            sort_by_distance(&mut asa_values, position);
            for asa_value in &asa_values {
                if asa_value.residue_id == residue.id().0 as i32 {
                    continue;
                }
                if asa_value.asa > exposed_asa_threshold {
                    let dpx = distance(&position, &asa_value.position);
                    dpx_residue_ids.push(residue.id().0 as i32);
                    dpx_values.push(dpx);
                    if args.pdb_output_path.is_some() {
                        for atom in residue.atoms_mut() {
                            atom.set_b_factor(dpx).unwrap();
                        }
                    }
                    break;
                }
            }
        } else {
            // Atom is exposed
            dpx_residue_ids.push(residue.id().0 as i32);
            dpx_values.push(0.0);
            if args.pdb_output_path.is_some() {
                for atom in residue.atoms_mut() {
                    atom.set_b_factor(0.0).unwrap();
                }
            }
        }
    }
    if let Some(pdb_path) = args.pdb_output_path {
        pdbtbx::save(&pdb, &pdb_path, pdbtbx::StrictnessLevel::Loose).unwrap();
        println!(
            "{}",
            format!(
                "Wrote output PDB to {} - DPX Values are saved in B-factor column",
                pdb_path
            )
            .green()
        );
    }
    if let Some(csv_path) = args.csv_output_path {
        let mut df: DataFrame = df!(
            "Residue ID" => &dpx_residue_ids,
            "DPX" => &dpx_values
        )
        .unwrap();
        let mut file = File::create(&csv_path).expect("could not create file");
        CsvWriter::new(&mut file)
            .include_header(true)
            .with_separator(b',')
            .finish(&mut df)
            .unwrap();
        println!("{} {}", "Wrote output CSV to ".green(), &csv_path);
    }
}
