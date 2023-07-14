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

use std::vec;

use crate::color_gradients::*;
use crate::pole_figures::{crystal_axis::CrystalAxes, minerals::Mineral};

use serde_derive::Deserialize;
/// The configuration of the pole figure analysis.
#[derive(Deserialize, Clone)]
pub struct PoleFiguresConfiguration {
    /// Optional value where to find the data file relating the numbers of the individual data files with the time they
    /// represent. The default value is `statistics`.
    #[serde(default = "time_data_file")]
    pub time_data_file: String,

    /// Optional value of the prefix where to find the data file containing global particle data per timestep such as the
    /// particle id, position, deformation type and elasticity information. The program will add a postfixes containing
    /// time and mpi process information in the format of `-00000.0000.dat` in which the first 5 zero's represent the timestep
    /// and the last 4 zero's represent the different files for that timestep. For each timestep, the program will automatically
    /// search through all files of the same timestep untill it find the file containing the particle id it was looking for.
    /// The program assumes that all the particles of the same id are not spread out over sever files, but stored in one
    /// file. Having multiple particle id's per file is fine.
    /// The default value is `particle_CPO/particles`.
    #[serde(default = "particle_data_file_prefix")]
    pub particle_data_file_prefix: String,

    /// Optional value of the prefix where to find the data file containing grain data per timestep. This is the particle
    /// id and euler angles of the olivine and enstatite. The program will add a postfixes containing
    /// time and mpi process information in the format of `-00000.0000.dat` in which the first 5 zero's represent the timestep
    /// and the last 4 zero's represent the different files for that timestep. For each timestep, the program will automatically
    /// search through all files of the same timestep untill it find the file containing the particle id it was looking for.
    /// The program assumes that all the particles of the same id are not spread out over sever files, but stored in one
    /// file. Having multiple particle id's per file is fine.
    /// The default value is `particle_CPO/weighted_CPO`.
    #[serde(default = "grain_data_file_prefix")]
    pub grain_data_file_prefix: String,

    /// Optional value of the prefix where to write out the produced pole figures. The program will add a postfix containing
    /// informatio about individual variables and settings such as the timestep and axes, etc.
    /// The default value is `CPO_figures/`.
    #[serde(default = "figure_output_dir")]
    pub figure_output_dir: String,

    /// Optional value of the prefix where to write out the produced pole figures. The program will add a postfix containing
    /// informatio about individual variables and settings such as the timestep and axes, etc.
    /// The default value is `weighted_CPO`.
    #[serde(default = "figure_output_prefix")]
    pub figure_output_prefix: String,

    /// Optional value of the color gradient. Available options are `Vik`, `Batlow`, `Simple`, `Imola`, `Hawaii` and `Roma`,
    /// which are all, except for `Simple`, from [http://www.fabiocrameri.ch/colourmaps.php](http://www.fabiocrameri.ch/colourmaps.php).
    /// The default value is `Batlow`.
    #[serde(default = "color_scale")]
    pub color_scale: ColorGradient,

    /// Optional value whether to inluce elasticity information in the header of the polefigure.
    /// Default is true.
    #[serde(default = "default_true")]
    pub elastisity_header: bool,

    /// Optional value whether to create a small (500x500 per pole figure) or normal (800x800 per pole figure) figure size.
    /// Default is false.
    #[serde(default = "default_false")]
    pub small_figure: bool,

    /// Optional value whether to omit mineral and axis information to the figure.
    /// Default value is false.
    #[serde(default = "default_false")]
    pub no_description_text: bool,

    /// A vector containing the times at which to make the pole figures. The acutal times are set to the closest time
    /// for which data is available.
    pub times: Vec<f64>,

    /// A vector containing the id's of the particles to plot. A plot for is made for every id.
    pub particle_ids: Vec<usize>,

    /// A vector of the crytal axis to plot. This will be plot on the horizontal axis of the plot, and the maximum
    /// of the maximum count of the pole figures will be used as the new maximum to scale the colors. Available options
    /// are `AAxis`, `BAxis` and `CAxis`.
    pub axes: Vec<CrystalAxes>,

    /// A vector containing the minerals to be plot. These will be added as a vertical axis to the plot. Available
    /// options are `Olivine` and `Enstatite`.
    pub minerals: Vec<Mineral>,
}

impl Default for PoleFiguresConfiguration {
    fn default() -> Self {
        PoleFiguresConfiguration {
            time_data_file: "statistics".to_string(),
            particle_data_file_prefix: "particles_cpo/particles".to_string(),
            grain_data_file_prefix: "particles_cpo/weighted_CPO".to_string(),
            figure_output_dir: "CPO_figures/".to_string(),
            figure_output_prefix: "weighted_CPO".to_string(),
            color_scale: ColorGradient::Batlow,
            elastisity_header: true,
            small_figure: false,
            no_description_text: false,
            times: vec![],
            particle_ids: vec![],
            axes: vec![],
            minerals: vec![],
        }
    }
}

fn default_false() -> bool {
    false
}

fn default_true() -> bool {
    true
}

fn time_data_file() -> String {
    PoleFiguresConfiguration {
        ..PoleFiguresConfiguration::default()
    }
    .time_data_file
}

fn particle_data_file_prefix() -> String {
    PoleFiguresConfiguration {
        ..PoleFiguresConfiguration::default()
    }
    .particle_data_file_prefix
}

fn grain_data_file_prefix() -> String {
    PoleFiguresConfiguration {
        ..PoleFiguresConfiguration::default()
    }
    .grain_data_file_prefix
}

fn figure_output_dir() -> String {
    PoleFiguresConfiguration {
        ..PoleFiguresConfiguration::default()
    }
    .figure_output_dir
}

fn figure_output_prefix() -> String {
    PoleFiguresConfiguration {
        ..PoleFiguresConfiguration::default()
    }
    .figure_output_prefix
}

fn color_scale() -> ColorGradient {
    PoleFiguresConfiguration {
        ..PoleFiguresConfiguration::default()
    }
    .color_scale
}
