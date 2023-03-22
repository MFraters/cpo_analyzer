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

use serde_derive::Deserialize;
/// A structure to load the CPO data
#[derive(Debug, Deserialize)]
#[allow(non_snake_case)]
pub struct Record {
    pub id: usize,
    pub mineral_0_EA_phi: Option<f64>,
    pub mineral_0_EA_theta: Option<f64>,
    pub mineral_0_EA_z: Option<f64>,
    pub mineral_1_EA_phi: Option<f64>,
    pub mineral_1_EA_theta: Option<f64>,
    pub mineral_1_EA_z: Option<f64>,
    pub mineral_0_dis_dens_0: Option<f64>,
    pub mineral_0_dis_dens_1: Option<f64>,
    pub mineral_0_dis_dens_2: Option<f64>,
    pub mineral_0_dis_dens_3: Option<f64>,
    pub mineral_0_recryst_frac_0: Option<f64>,
    pub mineral_0_recryst_frac_1: Option<f64>,
    pub mineral_0_recryst_frac_2: Option<f64>,
    pub mineral_0_recryst_frac_3: Option<f64>,
}
