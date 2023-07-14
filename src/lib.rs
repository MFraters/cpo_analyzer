/*
  Copyright (C) 2021 by the authors of the CPO Analyzer code.

  This file is part of the CPO Analyzer.

  The CPO Analyzer is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  The CPO Analyzer is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the CPO Analyzer; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

//! # CPO analyzer
//!
//! This crate contains tools to analyzer Crystal/Lattice Preffered Orientation (CPO/CPO) data. It is currently
//! only able to crate sets of polefigures from ASPECT CPO data, but other type of inputs and plots of analysis
//! are in the scope of this crate.
//!
//! Note that this crate is a very early (beta) release, and it functions for the specific purpose it was build
//! for, but a lot of functionality that is currently hard-coded can be generalized if there is interest. This
//! also means that the current interface and input file structure is subject to change.
//!
//! To run the analyzer, a configuration file is needed. The configuration files are written in the `.toml`
//! language. Here is an example file (where the lines starting with a `#` are comments):
//!
//! ```toml
//! # the location of the experiment dirs.
//! base_dir = "/path/to/base/dir/"
//!
//! # the directories containing the experiments. Currently only ASPECT output directories
//! # are supported.
//! experiment_dirs = ["experiment_1","experiment2"]
//!  
//! # Whether the Data was compressed with ZLIB.
//! compressed = true
//!
//! [pole_figures]
//!   # Wheter to include elasticity information in the header of the polefigure plots.
//!   elastisity_header = false
//!
//!   # For each time in this vector a new polefigure plot is made.
//!   times = [1.0,5.0,10]
//!
//!   # For each id in this vector a new polefigure plot is made.
//!   particle_ids = [1,10]
//!
//!   # A vector containing the pole figure axis to be plotted. These will be added as a
//!   # horizontal axis to the plot. Available options are `AAxis`, `BAxis` and `CAxis`.
//!   axes = ["AAxis","BAxis","CAxis"]
//!
//!   # A vector containing the minerals to be plotted. These will be added as a vertical
//!   # axis  to the plot. Available options are `Olivine` and `Enstatite`.
//!   minerals = ["Olivine","Enstatite"]
//! ```
//!
//! The configuration file without comments:
//! ```toml
//! base_dir = "/path/to/base/dir/"
//! experiment_dirs = ["experiment_1","experiment2"]
//! compressed = true
//!
//! [pole_figures]
//!   elastisity_header = false
//!   times = [1.0,5.0,10]
//!   particle_ids = [1,10]
//!   axes = ["AAxis","BAxis","CAxis"]
//!   minerals = ["Olivine","Enstatite"]
//! ```
pub mod color_gradients;
pub mod configuration;
pub mod pole_figures;

use crate::configuration::{
    config::Config, opt::Opt, particle_record::ParticleRecord, record::Record,
};
use crate::pole_figures::make_pole_figures::*;
use crate::pole_figures::{
    crystal_axis::CrystalAxes, lambert::*, minerals::Mineral, pole_figure::PoleFigure,
};

use ndarray::{Array, Array2, Axis};
use rayon::prelude::*;
use structopt::StructOpt;

use std::fs;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;

/// Entry point for the commandline agrument/binary version of the program. The location of the configuration file is taken
/// from those commandline arguments, and passed on to the `process` function.
///
/// # Example: CPO analyzer main function
///
/// ```should_panic
/// use cpo_analyzer::run;
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///   run()
/// }
/// ```
pub fn run() -> Result<(), Box<dyn std::error::Error>> {
    let opt = Opt::from_args();
    let config_file = opt.config_file;

    let configuration =
        load_configuration_file(config_file).expect("Could not load configuration file");
    process_configuration(configuration)
}

/// Load configuration TOML file into an internal structure which can be used by the `process_configuration` function.
///
/// # Example 1: CPO analyzer run function
///
/// ```
/// use structopt::StructOpt;
/// use cpo_analyzer::configuration::opt::Opt;
/// use cpo_analyzer::{load_configuration_file,process_configuration};
///
/// pub fn run() -> Result<(), Box<dyn std::error::Error>> {
///    let opt = Opt::from_args();
///    let config_file = opt.config_file;
///
///    let configuration = load_configuration_file(config_file).expect("Could not load configuration file");
///    process_configuration(configuration)
/// }
/// ```
///
/// # Example 2: stand-alone usage
///
/// ```
/// use std::path::PathBuf;
/// use cpo_analyzer::configuration::opt::Opt;
/// use cpo_analyzer::{load_configuration_file,process_configuration};
///
/// fn run() -> Result<(), Box<dyn std::error::Error>> {
///    let config_file = PathBuf::from("examples/config_example.toml");
///
///    let configuration = load_configuration_file(config_file).expect("Could not load configuration file");
///    process_configuration(configuration)
///}
/// ```
pub fn load_configuration_file(config_file: PathBuf) -> Result<Config, toml::de::Error> {
    let config_file_display = config_file.display();
    // Open the path in read-only mode, returns `io::Result<File>`
    let mut file = match File::open(&config_file) {
        // The `description` method of `io::Error` returns a string that
        // describes the error
        Err(why) => panic!("couldn't open {}: {}", config_file_display, why.to_string()),
        Ok(file) => file,
    };
    // Read the file contents into a string, returns `io::Result<usize>`
    let mut config_file_string = String::new();
    match file.read_to_string(&mut config_file_string) {
        Err(why) => panic!("couldn't read {}: {}", config_file_display, why.to_string()),
        Ok(_) => (),
    }

    toml::from_str(&config_file_string)
}

/// Get a list of time step from the statisitics file
fn get_time_steps(statistics_file: &String) -> Vec<f64> {
    println!("time data file:{}", statistics_file);
    let file = File::open(statistics_file).unwrap();
    let reader = BufReader::new(file);

    let mut data: String = "".to_string();
    for line in reader.lines() {
        let line = line.unwrap();
        let line = line.trim();
        let mut line = line.replace("  ", " ");
        while let Some(_) = line.find("  ") {
            line = line.replace("  ", " ");
        }

        if line.find('#') != Some(0) {
            if line.find("particles_cpo") != None {
                data = data + &line + "\n";
            }
        }
    }

    let mut time_step_to_time: Vec<f64> = vec![];
    let mut rdr = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .delimiter(b' ')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(data.as_bytes());
    for result in rdr.records() {
        // The iterator yields Result<StringRecord, Error>, so we check the
        // error here..
        let record = result.unwrap().clone();
        let time = record.get(1).clone();
        match time {
            Some(time) => time_step_to_time.push(time.to_string().parse::<f64>().unwrap()),
            None => assert!(false, "Time not found"),
        }
    }

    assert!(time_step_to_time.len() > 0,"Did not find any entries in the statistics file. This means there where no (non-comment) lines which had 'particle_CPO' in them.");

    time_step_to_time
}

/// find closest value to time in time_step_to_time
fn get_closest_time_step(output_time: &f64, time_step_to_time: &Vec<f64>) -> u64 {
    let after_time = time_step_to_time.iter().position(|x| x > &output_time);

    let after_time_step = match after_time {
        Some(time_step) => time_step,
        None => time_step_to_time.len() - 1,
    };

    // todo: oneline
    let before_time_step = if after_time_step > 0 {
        after_time_step - 1
    } else {
        after_time_step
    };

    // check wheter before_time_step or after_time_step is closer to output_time,
    // then use that one.
    let before_time_step_diff = (output_time - time_step_to_time[before_time_step]).abs();
    let after_time_step_diff = (output_time - time_step_to_time[after_time_step]).abs();

    let time_step = if before_time_step_diff < after_time_step_diff {
        before_time_step as u64
    } else {
        after_time_step as u64
    };

    time_step
}

/// Entry point for if the location of the config file is already know, such as the `run` function.
///
/// # Example 1: CPO analyzer run function
///
/// ```
/// use structopt::StructOpt;
/// use cpo_analyzer::configuration::opt::Opt;
/// use cpo_analyzer::{load_configuration_file,process_configuration};
///
/// pub fn run() -> Result<(), Box<dyn std::error::Error>> {
///    let opt = Opt::from_args();
///    let config_file = opt.config_file;
///
///    let configuration = load_configuration_file(config_file).expect("Could not load configuration file");
///    process_configuration(configuration)
/// }
/// ```
///
/// # Example 2: stand-alone usage
///
/// ```
/// use std::path::PathBuf;
/// use cpo_analyzer::configuration::opt::Opt;
/// use cpo_analyzer::{load_configuration_file,process_configuration};
///
/// fn run() -> Result<(), Box<dyn std::error::Error>> {
///    let config_file = PathBuf::from("examples/config_example.toml");
///
///    let configuration = load_configuration_file(config_file).expect("Could not load configuration file");
///    process_configuration(configuration)
///}
/// ```
pub fn process_configuration(config: Config) -> Result<(), Box<dyn std::error::Error>> {
    let before = Instant::now();

    let base_dir = config.base_dir.clone();

    // start the experimentsc
    let experiment_dirs = config.experiment_dirs.clone();
    experiment_dirs.par_iter().for_each(|experiment_dir| {
        println!("Processing experiment {}", experiment_dir);

        let cpo_dir = base_dir.clone() + &experiment_dir;

        if config.pole_figures.is_some() {
            // get a vector with the time for all the time steps
            let statistics_file =
                cpo_dir.to_owned() + &config.pole_figures.as_ref().unwrap().time_data_file;
            let time_step_to_time = get_time_steps(&statistics_file);
            let pole_figure_configuration = config.pole_figures.as_ref().unwrap();
            let elastisity_header = pole_figure_configuration.elastisity_header;

            for output_time in &pole_figure_configuration.times {
                // find closest value in time_step_to_time
                // assume it always starts a zero
                let time_step = get_closest_time_step(&output_time, &time_step_to_time);
                let time = time_step_to_time[time_step as usize];

                println!(
                    "Processing time {} (requested time: {}), located in time step : {}",
                    time, output_time, time_step,
                );

                fs::create_dir_all(
                    cpo_dir.to_owned() + &pole_figure_configuration.figure_output_dir,
                )
                .unwrap();

                let file_prefix_figures = pole_figure_configuration.figure_output_dir.to_owned()
                    + &pole_figure_configuration.figure_output_prefix;
                let mut rank_id = 0;

                let gam = 1.0; //0.5; // exponent for power-law normalization of color-scale
                               //let f = 1.05; // factor to make plot limits slightly bigger than the circle

                let sphere_points = 301; //76;//151;

                println!(
                    "particle ids size {}",
                    pole_figure_configuration.particle_ids.len()
                );
                for particle_id in &pole_figure_configuration.particle_ids {
                    println!("processing particle_id {}", particle_id);
                    let mut particle_olivine_a_axis_vectors = Vec::new();
                    let mut particle_olivine_b_axis_vectors = Vec::new();
                    let mut particle_olivine_c_axis_vectors = Vec::new();
                    let mut particle_enstatite_a_axis_vectors = Vec::new();
                    let mut particle_enstatite_b_axis_vectors = Vec::new();
                    let mut particle_enstatite_c_axis_vectors = Vec::new();

                    let mut file_found: bool = false;
                    while !file_found {
                        let angles_file = format!(
                            "{}{}-{:05}.{:04}.dat",
                            cpo_dir,
                            pole_figure_configuration.grain_data_file_prefix,
                            time_step,
                            rank_id
                        );
                        let angles_file = Path::new(&angles_file);

                        let mut config_mineral_string = String::new();
                        for mineral in pole_figure_configuration.minerals.clone() {
                            config_mineral_string = format!(
                                "{}{}",
                                config_mineral_string.clone(),
                                match mineral {
                                    Mineral::Olivine => {
                                        "oli_"
                                    }
                                    Mineral::Enstatite => {
                                        "ens_"
                                    }
                                }
                            )
                        }
                        let mut config_axis_string = String::new();
                        for axis in pole_figure_configuration.axes.clone() {
                            config_axis_string = format!(
                                "{}{}",
                                config_axis_string.clone(),
                                match axis {
                                    CrystalAxes::AAxis => {
                                        "A-"
                                    }
                                    CrystalAxes::BAxis => {
                                        "B-"
                                    }
                                    CrystalAxes::CAxis => {
                                        "C-"
                                    }
                                }
                            )
                        }
                        config_axis_string = format!("{}Axis_", config_axis_string.clone());

                        let output_file = format!(
                            "{}{}_{}{}{}{}_g{}_sp{}_t{:05}.{:05}.png",
                            cpo_dir,
                            file_prefix_figures,
                            if elastisity_header {
                                "elastic_"
                            } else {
                                "no-elastic_"
                            },
                            config_mineral_string,
                            config_axis_string,
                            pole_figure_configuration.color_scale.to_string(),
                            gam,
                            sphere_points,
                            time_step,
                            particle_id
                        );
                        let output_file = Path::new(&output_file);
                        let particle_file = format!(
                            "{}{}-{:05}.{:04}.dat",
                            cpo_dir,
                            pole_figure_configuration.particle_data_file_prefix,
                            time_step,
                            rank_id
                        );
                        let particle_info_file = Path::new(&particle_file);

                        println!("  trying file name: {}", angles_file.display());

                        // check wheter file exists, if not it means that is reached the max rank, so stop.
                        if !(fs::metadata(angles_file).is_ok()) {
                            println!(
                                "particle id {} not found for time step {}.",
                                particle_id, time_step
                            );
                            break;
                        }

                        // check wheter file is empty, if not continue to next rank
                        if fs::metadata(angles_file).unwrap().len() == 0 {
                            rank_id = rank_id + 1;
                            continue;
                        }

                        println!(
                            "  found particle id {} in:{}",
                            particle_id,
                            angles_file.display()
                        );
                        let file = File::open(angles_file).unwrap();
                        let metadata = file.metadata().unwrap();

                        let mut buf_reader =
                            BufReader::with_capacity(metadata.len() as usize, file);

                        let mut decoded_data = Vec::new();

                        let compressed = config.compressed;

                        let decoded_reader = if compressed {
                            let mut decoder = libflate::zlib::Decoder::new(buf_reader).unwrap();
                            decoder.read_to_end(&mut decoded_data).unwrap();
                            String::from_utf8_lossy(&decoded_data)
                        } else {
                            let data = buf_reader.fill_buf().unwrap();
                            String::from_utf8_lossy(&data)
                        };

                        let mut rdr = csv::ReaderBuilder::new()
                            .has_headers(true)
                            .delimiter(b' ')
                            .from_reader(decoded_reader.as_bytes());

                        let mut integer = 0;
                        for result in rdr.deserialize() {
                            println!("result: {:?}", result);
                            let record: Record = result.unwrap();
                            if record.id == *particle_id {
                                let deg_to_rad = std::f64::consts::PI / 180.;

                                // olivine
                                if pole_figure_configuration
                                    .minerals
                                    .contains(&Mineral::Olivine)
                                {
                                    let euler_angles = Array::from(vec![
                                        record.mineral_0_EA_phi.unwrap() * deg_to_rad,
                                        record.mineral_0_EA_theta.unwrap() * deg_to_rad,
                                        record.mineral_0_EA_z.unwrap() * deg_to_rad,
                                    ]);
                                    let rotation_matrix =
                                        euler_angles_to_rotation_matrix(euler_angles).unwrap();

                                    particle_olivine_a_axis_vectors
                                        .push(rotation_matrix.row(0).to_owned());
                                    particle_olivine_b_axis_vectors
                                        .push(rotation_matrix.row(1).to_owned());
                                    particle_olivine_c_axis_vectors
                                        .push(rotation_matrix.row(2).to_owned());
                                }

                                // enstatite
                                if pole_figure_configuration
                                    .minerals
                                    .contains(&Mineral::Enstatite)
                                {
                                    let euler_angles = Array::from(vec![
                                        record.mineral_1_EA_phi.unwrap() * deg_to_rad,
                                        record.mineral_1_EA_theta.unwrap() * deg_to_rad,
                                        record.mineral_1_EA_z.unwrap() * deg_to_rad,
                                    ]);
                                    let rotation_matrix =
                                        euler_angles_to_rotation_matrix(euler_angles).unwrap();

                                    particle_enstatite_a_axis_vectors
                                        .push(rotation_matrix.row(0).to_owned());
                                    particle_enstatite_b_axis_vectors
                                        .push(rotation_matrix.row(1).to_owned());
                                    particle_enstatite_c_axis_vectors
                                        .push(rotation_matrix.row(2).to_owned());
                                }
                            }
                            integer = integer + 1;
                        }

                        // check if the particle id was found in this file, otherwise continue
                        if particle_olivine_a_axis_vectors.len() == 0
                            && particle_enstatite_a_axis_vectors.len() == 0
                        {
                            rank_id = rank_id + 1;
                            continue;
                        }
                        file_found = true;

                        // retrieve anisotropy info
                        let mut particle_record = ParticleRecord {
                            id: 0,
                            x: 0.0,
                            y: 0.0,
                            z: Some(0.0),
                            olivine_deformation_type: None,
                            full_norm_square: None,
                            triclinic_norm_square_p1: None,
                            triclinic_norm_square_p2: None,
                            triclinic_norm_square_p3: None,
                            monoclinic_norm_square_p1: None,
                            monoclinic_norm_square_p2: None,
                            monoclinic_norm_square_p3: None,
                            orthohombic_norm_square_p1: None,
                            orthohombic_norm_square_p2: None,
                            orthohombic_norm_square_p3: None,
                            tetragonal_norm_square_p1: None,
                            tetragonal_norm_square_p2: None,
                            tetragonal_norm_square_p3: None,
                            hexagonal_norm_square_p1: None,
                            hexagonal_norm_square_p2: None,
                            hexagonal_norm_square_p3: None,
                            isotropic_norm_square: None,
                        };

                        let particle_info_file = File::open(particle_info_file).unwrap();
                        let buf_reader = BufReader::new(particle_info_file);

                        let mut rdr = csv::ReaderBuilder::new()
                            .has_headers(true)
                            .delimiter(b' ')
                            .from_reader(buf_reader);

                        for result in rdr.deserialize() {
                            // We must tell Serde what type we want to deserialize into.
                            let record: ParticleRecord = result.unwrap();
                            if record.id == *particle_id {
                                particle_record = record;
                            }
                        }
                        // end retrieve anisotropy info
                        println!("end retrieve antisotropy info");
                        println!("create lambert equal area gridpoint");

                        let lambert =
                            create_lambert_equal_area_gridpoint(sphere_points, "upper".to_string())
                                .unwrap();

                        println!("create sphere_point_grid");
                        let mut sphere_point_grid =
                            Array2::zeros((3, sphere_points * sphere_points));

                        for i in 0..sphere_points {
                            for j in 0..sphere_points {
                                sphere_point_grid[[0, i * sphere_points + j]] = lambert.x[[i, j]];
                                sphere_point_grid[[1, i * sphere_points + j]] = lambert.y[[i, j]];
                                sphere_point_grid[[2, i * sphere_points + j]] = lambert.z[[i, j]];
                            }
                        }

                        let n_grains = std::cmp::max(
                            particle_olivine_a_axis_vectors.len(),
                            particle_enstatite_a_axis_vectors.len(),
                        );

                        let mut pole_figure_grid: Vec<Vec<PoleFigure>> = vec![
                                vec![
                                    PoleFigure {
                                        crystal_axis: CrystalAxes::AAxis,
                                        mineral: Mineral::Olivine,
                                        counts: Array2::zeros((n_grains, 3)),
                                        max_count: 0.0,
                                    };
                                    pole_figure_configuration.minerals.len()
                                ];
                                pole_figure_configuration.axes.len()
                            ];

                        let mut figure_horizontal_axis = 0;
                        for axis in pole_figure_configuration.axes.clone() {
                            let mut figure_vertical_axis = 0;
                            for mineral in pole_figure_configuration.minerals.clone() {
                                let mut particle_arrays = Array2::zeros((n_grains, 3));
                                for i in 0..n_grains {
                                    for j in 0..3 {
                                        particle_arrays[[i, j]] = match axis {
                                            CrystalAxes::AAxis => match mineral {
                                                Mineral::Olivine => {
                                                    particle_olivine_a_axis_vectors[i][j]
                                                }
                                                Mineral::Enstatite => {
                                                    particle_enstatite_a_axis_vectors[i][j]
                                                }
                                            },
                                            CrystalAxes::BAxis => match mineral {
                                                Mineral::Olivine => {
                                                    particle_olivine_b_axis_vectors[i][j]
                                                }
                                                Mineral::Enstatite => {
                                                    particle_enstatite_b_axis_vectors[i][j]
                                                }
                                            },
                                            CrystalAxes::CAxis => match mineral {
                                                Mineral::Olivine => {
                                                    particle_olivine_c_axis_vectors[i][j]
                                                }
                                                Mineral::Enstatite => {
                                                    particle_enstatite_c_axis_vectors[i][j]
                                                }
                                            },
                                        };
                                    }
                                }
                                let counts = gaussian_orientation_counts(
                                    &particle_arrays,
                                    &sphere_point_grid,
                                    sphere_points,
                                )
                                .unwrap();

                                let mut max_count_value = 0.0;

                                for i in 0..counts.shape()[0] - 1 {
                                    for j in 0..counts.shape()[1] - 1 {
                                        if counts[[i, j]] > max_count_value {
                                            max_count_value = counts[[i, j]];
                                        }
                                    }
                                }

                                pole_figure_grid[figure_horizontal_axis][figure_vertical_axis] =
                                    PoleFigure {
                                        crystal_axis: axis.clone(),
                                        mineral: mineral.clone(),
                                        counts: counts,
                                        max_count: max_count_value,
                                    };
                                figure_vertical_axis += 1;
                            }
                            figure_horizontal_axis += 1;
                        }

                        // set all horizontal max values to the max of the horizontal max max values
                        let mut figure_vertical_axis = 0;
                        for _mineral in pole_figure_configuration.minerals.clone() {
                            let mut max_count_value = 0.0;

                            // loop to find the max value
                            let mut figure_horizontal_axis = 0;
                            for _axis in pole_figure_configuration.axes.clone() {
                                if pole_figure_grid[figure_horizontal_axis][figure_vertical_axis]
                                    .max_count
                                    > max_count_value
                                {
                                    max_count_value = pole_figure_grid[figure_horizontal_axis]
                                        [figure_vertical_axis]
                                        .max_count;
                                }
                                figure_horizontal_axis += 1;
                            }

                            // loop to write max value
                            figure_horizontal_axis = 0;
                            for _axis in pole_figure_configuration.axes.clone() {
                                pole_figure_grid[figure_horizontal_axis][figure_vertical_axis]
                                    .max_count = max_count_value;
                                figure_horizontal_axis += 1;
                            }
                            figure_vertical_axis += 1;
                        }

                        make_pole_figures(
                            pole_figure_configuration.small_figure,
                            pole_figure_configuration.no_description_text,
                            elastisity_header,
                            n_grains,
                            0,
                            &pole_figure_grid,
                            &lambert,
                            output_file,
                            &particle_record,
                            time,
                            gam,
                            &pole_figure_configuration.color_scale,
                        )
                        .unwrap();

                        println!(
                            "  After make_polefigures: Elapsed time: {:.2?}",
                            before.elapsed()
                        );
                    }
                    println!("go to next id");
                }
            }
        }
    });
    Ok(())
}

/// Utility function to compute a rotation matrix from Z-X-Z Euler angles.
fn euler_angles_to_rotation_matrix(
    euler_angles: Array<f64, ndarray::Dim<[usize; 1]>>, //phi1: f64,
                                                        //theta: f64,
                                                        //phi2: f64
) -> Result<Array2<f64>, Box<dyn std::error::Error>> {
    let mut rotation_matrix: Array2<f64> = Array::zeros((3, 3));

    rotation_matrix[[0, 0]] = euler_angles[2].cos() * euler_angles[0].cos()
        - euler_angles[1].cos() * euler_angles[0].sin() * euler_angles[2].sin();
    rotation_matrix[[0, 1]] = -euler_angles[2].cos() * euler_angles[0].sin()
        - euler_angles[1].cos() * euler_angles[0].cos() * euler_angles[2].sin();
    rotation_matrix[[0, 2]] = -euler_angles[2].sin() * euler_angles[1].sin();

    rotation_matrix[[1, 0]] = euler_angles[2].sin() * euler_angles[0].cos()
        + euler_angles[1].cos() * euler_angles[0].sin() * euler_angles[2].cos();
    rotation_matrix[[1, 1]] = -euler_angles[2].sin() * euler_angles[0].sin()
        + euler_angles[1].cos() * euler_angles[0].cos() * euler_angles[2].cos();
    rotation_matrix[[1, 2]] = euler_angles[2].cos() * euler_angles[1].sin();

    rotation_matrix[[2, 0]] = -euler_angles[1].sin() * euler_angles[0].sin();
    rotation_matrix[[2, 1]] = -euler_angles[1].sin() * euler_angles[0].cos();
    rotation_matrix[[2, 2]] = euler_angles[1].cos();

    Ok(rotation_matrix)
}

/// Following method in Robin and Jowett, Tectonophysics, 1986
/// Computerized contouring and statistical evaluation of orientation data
/// using contouring circles and continuous weighting functions.
/// For the k value we use a combination between option 2 and 3, where option
/// 2 is used as long as k is larger as 100, otherwise it is set to 100.
fn gaussian_orientation_counts(
    particles: &Array2<f64>,
    sphere_point_grid: &Array2<f64>,
    sphere_points: usize,
) -> Result<Array2<f64>, Box<dyn std::error::Error>> {
    let npts = particles.shape()[0];

    // Choose k, which defines width of spherical gaussian  (table 3)
    let k = (2. * (1. + npts as f64 / 9.)).min(100.);

    // Given k, calculate standard deviation (eq 13b)
    let std_dev = ((npts as f64 * (k as f64 / 2. - 1.) / (k * k)) as f64).sqrt();

    // Calculate dot product
    let mut cosalpha = particles.dot(sphere_point_grid);

    // Calculate the counts from the spherical gaussian
    //let counts = Array::zeros(cosalpha.shape());
    cosalpha.par_mapv_inplace(f64::abs);

    cosalpha = (k as f64) * (cosalpha - 1.);

    cosalpha.par_mapv_inplace(f64::exp);

    let counts = cosalpha.sum_axis(Axis(0));
    let counts = counts.into_shape((sphere_points, sphere_points))?;

    // normalize so each MUD is 3 sigma from that expected for a uniform
    // distribution
    let counts = counts / (3. * std_dev);

    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn euler_angles_from_rotation_matrix(
        rotation_matrix: Array2<f64>,
    ) -> Array<f64, ndarray::Dim<[usize; 1]>> {
        // ZXZ Euler angles
        let mut euler_angles = Array::zeros(3);

        let theta = rotation_matrix[[2, 2]].acos();
        let mut phi1 = 0.0;
        let phi2;

        if theta != 0.0 && theta != std::f64::consts::PI {
            phi1 = (rotation_matrix[[2, 0]] / -theta.sin())
                .atan2(rotation_matrix[[2, 1]] / -theta.sin());
            phi2 = (rotation_matrix[[0, 2]] / -theta.sin())
                .atan2(rotation_matrix[[1, 2]] / theta.sin());
        } else {
            // note that in the case theta is 0 or phi a dimension is lost
            // see: https://en.wikipedia.org/wiki/Gimbal_lock. We set phi1
            // to 0 and compute the corresponding phi2. The resulting direction
            // (cosine matrix) should be the same.
            if theta == 0.0 {
                phi2 = -phi1 - rotation_matrix[[0, 1]].atan2(rotation_matrix[[0, 0]]);
            } else {
                phi2 = phi1 + rotation_matrix[[0, 1]].atan2(rotation_matrix[[0, 0]]);
            }
        }

        euler_angles[0] = phi1;
        euler_angles[1] = theta;
        euler_angles[2] = phi2;

        euler_angles[0] = euler_angles[0]
            - 2.0 * std::f64::consts::PI * (euler_angles[0] / (2.0 * std::f64::consts::PI)).floor();
        euler_angles[1] = euler_angles[1]
            - 2.0 * std::f64::consts::PI * (euler_angles[1] / (2.0 * std::f64::consts::PI)).floor();
        euler_angles[2] = euler_angles[2]
            - 2.0 * std::f64::consts::PI * (euler_angles[2] / (2.0 * std::f64::consts::PI)).floor();

        euler_angles
    }

    #[test]
    fn test_euler_angle_and_rotation_matrix_functions_part_1() {
        let mut rot1 = Array2::zeros((3, 3));
        rot1[[0, 0]] = 0.36;
        rot1[[0, 1]] = 0.48;
        rot1[[0, 2]] = -0.8;

        rot1[[1, 0]] = -0.8;
        rot1[[1, 1]] = 0.6;
        rot1[[1, 2]] = 0.;

        rot1[[2, 0]] = 0.48;
        rot1[[2, 1]] = 0.64;
        rot1[[2, 2]] = 0.6;

        let ea1 = euler_angles_from_rotation_matrix(rot1.clone());
        let rot2 = euler_angles_to_rotation_matrix(ea1).unwrap();

        assert!(
            (rot1[[0, 0]] - rot2[[0, 0]]).abs() < 1e-8,
            "Difference between rot1[0,0] ({}) and rot2[0,0] ({}) is too large: {}.",
            rot1[[0, 0]],
            rot2[[0, 0]],
            (rot1[[0, 0]] - rot2[[0, 0]]).abs()
        );
        assert!(
            (rot1[[0, 1]] - rot2[[0, 1]]).abs() < 1e-8,
            "Difference between rot1[0,1] ({}) and rot2[0,1] ({}) is too large: {}.",
            rot1[[0, 1]],
            rot2[[0, 1]],
            (rot1[[0, 1]] - rot2[[0, 1]]).abs()
        );
        assert!(
            (rot1[[0, 2]] - rot2[[0, 2]]).abs() < 1e-8,
            "Difference between rot1[0,2] ({}) and rot2[0,2] ({}) is too large: {}.",
            rot1[[0, 2]],
            rot2[[0, 2]],
            (rot1[[0, 2]] - rot2[[0, 2]]).abs()
        );
        assert!(
            (rot1[[1, 0]] - rot2[[1, 0]]).abs() < 1e-8,
            "Difference between rot1[1,0] ({}) and rot2[1,0] ({}) is too large: {}.",
            rot1[[1, 0]],
            rot2[[1, 0]],
            (rot1[[1, 0]] - rot2[[1, 0]]).abs()
        );
        assert!(
            (rot1[[1, 1]] - rot2[[1, 1]]).abs() < 1e-8,
            "Difference between rot1[1,1] ({}) and rot2[1,1] ({}) is too large: {}.",
            rot1[[1, 1]],
            rot2[[1, 1]],
            (rot1[[1, 1]] - rot2[[1, 1]]).abs()
        );
        assert!(
            (rot1[[1, 2]] - rot2[[1, 2]]).abs() < 1e-8,
            "Difference between rot1[1,2] ({}) and rot2[1,2] ({}) is too large: {}.",
            rot1[[1, 2]],
            rot2[[1, 2]],
            (rot1[[1, 2]] - rot2[[1, 2]]).abs()
        );
        assert!(
            (rot1[[2, 0]] - rot2[[2, 0]]).abs() < 1e-8,
            "Difference between rot1[2,0] ({}) and rot2[2,0] ({}) is too large: {}.",
            rot1[[2, 0]],
            rot2[[2, 0]],
            (rot1[[2, 0]] - rot2[[2, 0]]).abs()
        );
        assert!(
            (rot1[[2, 1]] - rot2[[2, 1]]).abs() < 1e-8,
            "Difference between rot1[2,1] ({}) and rot2[2,1] ({}) is too large: {}.",
            rot1[[2, 1]],
            rot2[[2, 1]],
            (rot1[[2, 1]] - rot2[[2, 1]]).abs()
        );
        assert!(
            (rot1[[2, 2]] - rot2[[2, 2]]).abs() < 1e-8,
            "Difference between rot1[2,2] ({}) and rot2[2,2] ({}) is too large: {}.",
            rot1[[2, 2]],
            rot2[[2, 2]],
            (rot1[[2, 2]] - rot2[[2, 2]]).abs()
        );

        let ea2 = euler_angles_from_rotation_matrix(rot2.clone());
        let rot3 = euler_angles_to_rotation_matrix(ea2).unwrap();

        assert!(
            (rot1[[0, 0]] - rot3[[0, 0]]).abs() < 1e-8,
            "Difference between rot1[0,0] ({}) and rot3[0,0] ({}) is too large: {}.",
            rot1[[0, 0]],
            rot3[[0, 0]],
            (rot1[[0, 0]] - rot3[[0, 0]]).abs()
        );
        assert!(
            (rot1[[0, 1]] - rot3[[0, 1]]).abs() < 1e-8,
            "Difference between rot1[0,1] ({}) and rot3[0,1] ({}) is too large: {}.",
            rot1[[0, 1]],
            rot3[[0, 1]],
            (rot1[[0, 1]] - rot3[[0, 1]]).abs()
        );
        assert!(
            (rot1[[0, 2]] - rot3[[0, 2]]).abs() < 1e-8,
            "Difference between rot1[0,2] ({}) and rot3[0,2] ({}) is too large: {}.",
            rot1[[0, 2]],
            rot3[[0, 2]],
            (rot1[[0, 2]] - rot3[[0, 2]]).abs()
        );
        assert!(
            (rot1[[1, 0]] - rot3[[1, 0]]).abs() < 1e-8,
            "Difference between rot1[1,0] ({}) and rot3[1,0] ({}) is too large: {}.",
            rot1[[1, 0]],
            rot3[[1, 0]],
            (rot1[[1, 0]] - rot3[[1, 0]]).abs()
        );
        assert!(
            (rot1[[1, 1]] - rot3[[1, 1]]).abs() < 1e-8,
            "Difference between rot1[1,1] ({}) and rot3[1,1] ({}) is too large: {}.",
            rot1[[1, 1]],
            rot3[[1, 1]],
            (rot1[[1, 1]] - rot3[[1, 1]]).abs()
        );
        assert!(
            (rot1[[1, 2]] - rot3[[1, 2]]).abs() < 1e-8,
            "Difference between rot1[1,2] ({}) and rot3[1,2] ({}) is too large: {}.",
            rot1[[1, 2]],
            rot3[[1, 2]],
            (rot1[[1, 2]] - rot3[[1, 2]]).abs()
        );
        assert!(
            (rot1[[2, 0]] - rot3[[2, 0]]).abs() < 1e-8,
            "Difference between rot1[2,0] ({}) and rot3[2,0] ({}) is too large: {}.",
            rot1[[2, 0]],
            rot3[[2, 0]],
            (rot1[[2, 0]] - rot3[[2, 0]]).abs()
        );
        assert!(
            (rot1[[2, 1]] - rot3[[2, 1]]).abs() < 1e-8,
            "Difference between rot1[2,1] ({}) and rot3[2,1] ({}) is too large: {}.",
            rot1[[2, 1]],
            rot3[[2, 1]],
            (rot1[[2, 1]] - rot3[[2, 1]]).abs()
        );
        assert!(
            (rot1[[2, 2]] - rot3[[2, 2]]).abs() < 1e-8,
            "Difference between rot1[2,2] ({}) and rot3[2,2] ({}) is too large: {}.",
            rot1[[2, 2]],
            rot3[[2, 2]],
            (rot1[[2, 2]] - rot3[[2, 2]]).abs()
        );
    }

    #[test]
    fn test_euler_angle_and_rotation_matrix_functions_part_2() {
        let mut rot1 = Array2::zeros((3, 3));
        rot1[[0, 0]] = 0.36;
        rot1[[0, 1]] = 0.48;
        rot1[[0, 2]] = -0.8;

        rot1[[1, 0]] = -0.8;
        rot1[[1, 1]] = 0.6;
        rot1[[1, 2]] = 0.;

        rot1[[2, 0]] = 0.0;
        rot1[[2, 1]] = 0.0;
        rot1[[2, 2]] = 0.0;

        let mut rot2_expected = Array2::zeros((3, 3));
        rot2_expected[[0, 0]] = 0.;
        rot2_expected[[0, 1]] = 0.;
        rot2_expected[[0, 2]] = -1.;
        rot2_expected[[1, 0]] = -1.;
        rot2_expected[[1, 1]] = 0.0;
        rot2_expected[[1, 2]] = 0.0;
        rot2_expected[[2, 0]] = 0.0;
        rot2_expected[[2, 1]] = 1.0;
        rot2_expected[[2, 2]] = 0.0;

        let ea1 = euler_angles_from_rotation_matrix(rot1.clone());
        let rot2 = euler_angles_to_rotation_matrix(ea1).unwrap();

        assert!(
            (rot2_expected[[0, 0]] - rot2[[0, 0]]).abs() < 1e-8,
            "Difference between rot2_expected[0,0] ({}) and rot2[0,0] ({}) is too large: {}.",
            rot2_expected[[0, 0]],
            rot2[[0, 0]],
            (rot2_expected[[0, 0]] - rot2[[0, 0]]).abs()
        );
        assert!(
            (rot2_expected[[0, 1]] - rot2[[0, 1]]).abs() < 1e-8,
            "Difference between rot2_expected[0,1] ({}) and rot2[0,1] ({}) is too large: {}.",
            rot2_expected[[0, 1]],
            rot2[[0, 1]],
            (rot2_expected[[0, 1]] - rot2[[0, 1]]).abs()
        );
        assert!(
            (rot2_expected[[0, 2]] - rot2[[0, 2]]).abs() < 1e-8,
            "Difference between rot2_expected[0,2] ({}) and rot2[0,2] ({}) is too large: {}.",
            rot2_expected[[0, 2]],
            rot2[[0, 2]],
            (rot2_expected[[0, 2]] - rot2[[0, 2]]).abs()
        );
        assert!(
            (rot2_expected[[1, 0]] - rot2[[1, 0]]).abs() < 1e-8,
            "Difference between rot2_expected[1,0] ({}) and rot2[1,0] ({}) is too large: {}.",
            rot2_expected[[1, 0]],
            rot2[[1, 0]],
            (rot2_expected[[1, 0]] - rot2[[1, 0]]).abs()
        );
        assert!(
            (rot2_expected[[1, 1]] - rot2[[1, 1]]).abs() < 1e-8,
            "Difference between rot2_expected[1,1] ({}) and rot2[1,1] ({}) is too large: {}.",
            rot2_expected[[1, 1]],
            rot2[[1, 1]],
            (rot2_expected[[1, 1]] - rot2[[1, 1]]).abs()
        );
        assert!(
            (rot2_expected[[1, 2]] - rot2[[1, 2]]).abs() < 1e-8,
            "Difference between rot2_expected[1,2] ({}) and rot2[1,2] ({}) is too large: {}.",
            rot2_expected[[1, 2]],
            rot2[[1, 2]],
            (rot2_expected[[1, 2]] - rot2[[1, 2]]).abs()
        );
        assert!(
            (rot2_expected[[2, 0]] - rot2[[2, 0]]).abs() < 1e-8,
            "Difference between rot2_expected[2,0] ({}) and rot2[2,0] ({}) is too large: {}.",
            rot2_expected[[2, 0]],
            rot2[[2, 0]],
            (rot2_expected[[2, 0]] - rot2[[2, 0]]).abs()
        );
        assert!(
            (rot2_expected[[2, 1]] - rot2[[2, 1]]).abs() < 1e-8,
            "Difference between rot2_expected[2,1] ({}) and rot2[2,1] ({}) is too large: {}.",
            rot2_expected[[2, 1]],
            rot2[[2, 1]],
            (rot2_expected[[2, 1]] - rot2[[2, 1]]).abs()
        );
        assert!(
            (rot2_expected[[2, 2]] - rot2[[2, 2]]).abs() < 1e-8,
            "Difference between rot2_expected[2,2] ({}) and rot2[2,2] ({}) is too large: {}.",
            rot2_expected[[2, 2]],
            rot2[[2, 2]],
            (rot2_expected[[2, 2]] - rot2[[2, 2]]).abs()
        );

        let ea2 = euler_angles_from_rotation_matrix(rot2.clone());
        let rot3 = euler_angles_to_rotation_matrix(ea2).unwrap();

        assert!(
            (rot2_expected[[0, 0]] - rot3[[0, 0]]).abs() < 1e-8,
            "Difference between rot2_expected[0,0] ({}) and rot3[0,0] ({}) is too large: {}.",
            rot2_expected[[0, 0]],
            rot3[[0, 0]],
            (rot2_expected[[0, 0]] - rot3[[0, 0]]).abs()
        );
        assert!(
            (rot2_expected[[0, 1]] - rot3[[0, 1]]).abs() < 1e-8,
            "Difference between rot2_expected[0,1] ({}) and rot3[0,1] ({}) is too large: {}.",
            rot2_expected[[0, 1]],
            rot3[[0, 1]],
            (rot2_expected[[0, 1]] - rot3[[0, 1]]).abs()
        );
        assert!(
            (rot2_expected[[0, 2]] - rot3[[0, 2]]).abs() < 1e-8,
            "Difference between rot2_expected[0,2] ({}) and rot3[0,2] ({}) is too large: {}.",
            rot2_expected[[0, 2]],
            rot3[[0, 2]],
            (rot2_expected[[0, 2]] - rot3[[0, 2]]).abs()
        );
        assert!(
            (rot2_expected[[1, 0]] - rot3[[1, 0]]).abs() < 1e-8,
            "Difference between rot2_expected[1,0] ({}) and rot3[1,0] ({}) is too large: {}.",
            rot2_expected[[1, 0]],
            rot3[[1, 0]],
            (rot2_expected[[1, 0]] - rot3[[1, 0]]).abs()
        );
        assert!(
            (rot2_expected[[1, 1]] - rot3[[1, 1]]).abs() < 1e-8,
            "Difference between rot2_expected[1,1] ({}) and rot3[1,1] ({}) is too large: {}.",
            rot2_expected[[1, 1]],
            rot3[[1, 1]],
            (rot2_expected[[1, 1]] - rot3[[1, 1]]).abs()
        );
        assert!(
            (rot2_expected[[1, 2]] - rot3[[1, 2]]).abs() < 1e-8,
            "Difference between rot2_expected[1,2] ({}) and rot3[1,2] ({}) is too large: {}.",
            rot2_expected[[1, 2]],
            rot3[[1, 2]],
            (rot2_expected[[1, 2]] - rot3[[1, 2]]).abs()
        );
        assert!(
            (rot2_expected[[2, 0]] - rot3[[2, 0]]).abs() < 1e-8,
            "Difference between rot2_expected[2,0] ({}) and rot3[2,0] ({}) is too large: {}.",
            rot2_expected[[2, 0]],
            rot3[[2, 0]],
            (rot2_expected[[2, 0]] - rot3[[2, 0]]).abs()
        );
        assert!(
            (rot2_expected[[2, 1]] - rot3[[2, 1]]).abs() < 1e-8,
            "Difference between rot2_expected[2,1] ({}) and rot3[2,1] ({}) is too large: {}.",
            rot2_expected[[2, 1]],
            rot3[[2, 1]],
            (rot2_expected[[2, 1]] - rot3[[2, 1]]).abs()
        );
        assert!(
            (rot2_expected[[2, 2]] - rot3[[2, 2]]).abs() < 1e-8,
            "Difference between rot2_expected[2,2] ({}) and rot3[2,2] ({}) is too large: {}.",
            rot2_expected[[2, 2]],
            rot3[[2, 2]],
            (rot2_expected[[2, 2]] - rot3[[2, 2]]).abs()
        );
    }
}
