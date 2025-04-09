#!/usr/bin/env python3
## Alireza Rashti - Apr 2025 (C)
## usage:
## $ ./me -h
##

import argparse

# import math as m
# import numpy as np
# import sys
# import os
# import re
# import matplotlib.pyplot as plt
# import glob
# import sympy
## ---------------------------------------------------------------------- ##

# dummy par file
g_dummy = """
#-----------------------------------------------------------------------#
# Physics:                                                              #
#-----------------------------------------------------------------------#

#### Project:
Project = NS_NS_binary_initial_data

#### NS-NS binary:
NSNS_separation                  = {NSNS_separation}
NSNS_angular_velocity            = auto
NSNS_start_off                   = parameter_file
NSNS_observe_ADM_P               = S_obj1+S_obj2,default
NSNS_observe_ADM_J               = S_obj1+S_obj2,default
NSNS_P_ADM_control_method        = adjust(x_CM,y_CM)
NSNS_P_ADM_control_update_weight = 0.(x10)->0.2
NSNS_P_ADM_control_tolerance     = 1E-8
NSNS_P_ADM_control_threshold     = 10


#### NS1:
NS1_baryonic_mass               = {NS1_baryonic_mass}
NS1_EoS_description             = custom_tab_eos
NS1_EoS_type                    = tabular
NS1_EoS_unit                    = geo
NS1_EoS_table_path              = {NS1_EoS_table_path}
NS1_EoS_table_format            = rest_mass_density,specific_internal_energy,pressure
NS1_EoS_interpolation_method    = Hermite1D
NS1_EoS_interpolation_use_log   = yes
NS1_EoS_Hermite1D_FD_accuracy   = 3
NS1_EoS_Hermite1D_num_points    = 2
NS1_EoS_enthalpy_floor          = {NS1_EoS_enthalpy_floor}
NS1_EoS_enthalpy_ceiling        = {NS1_EoS_enthalpy_ceiling}
NS1_Omega_x                     = 0.
NS1_Omega_y                     = 0.
NS1_Omega_z                     = 0.
NS1_surface_type                = perfect_s2->topology_s2
NS1_surface_finder              = bisection
NS1_surface_change_threshold    = 0.0
NS1_surface_Ylm_max_l           = 10
NS1_enthalpy_allowed_residual   = 1E-8
NS1_enthalpy_update_weight      = 0.5 ## aggressive update!
NS1_Euler_const_criteria        = fix_baryonic_mass
NS1_Euler_const_update_weight   = 1.
NS1_force_balance_equation      = none(x4)->adjust(d/dy:Omega)
NS1_force_balance_update_weight = 0.2
NS1_adjust_center_method        = taylor_expansion
NS1_adjust_center_update_weight = 1.
NS1_extrapolate_matter_fields   = inverse_r_expmr
NS1_Eq_phi_polish               = 0.1
NS1_start_off                   = TOV

#### NS2:
NS2_baryonic_mass               = {NS2_baryonic_mass}
NS2_EoS_description             = custom_tab_eos
NS2_EoS_type                    = tabular
NS2_EoS_unit                    = geo
NS2_EoS_table_path              = {NS2_EoS_table_path}
NS2_EoS_table_format            = rest_mass_density,specific_internal_energy,pressure
NS2_EoS_interpolation_method    = Hermite1D
NS2_EoS_interpolation_use_log   = yes
NS2_EoS_Hermite1D_FD_accuracy   = 3
NS2_EoS_Hermite1D_num_points    = 2
NS2_EoS_enthalpy_floor          = {NS2_EoS_enthalpy_floor}
NS2_EoS_enthalpy_ceiling        = {NS2_EoS_enthalpy_ceiling}
NS2_Omega_x                     = 0.
NS2_Omega_y                     = 0.
NS2_Omega_z                     = 0.
NS2_surface_type                = perfect_s2->topology_s2
NS2_surface_finder              = bisection
NS2_surface_change_threshold    = 0.0
NS2_surface_Ylm_max_l           = 10
NS2_enthalpy_allowed_residual   = 1E-8
NS2_enthalpy_update_weight      = 0.5 ## aggressive update!
NS2_Euler_const_criteria        = fix_baryonic_mass
NS2_Euler_const_update_weight   = 1.
NS2_force_balance_equation      = none(x4)->adjust(d/dy:Omega)
NS2_force_balance_update_weight = 0.2
NS2_adjust_center_method        = taylor_expansion
NS2_adjust_center_update_weight = 1.
NS2_extrapolate_matter_fields   = inverse_r_expmr
NS2_Eq_phi_polish               = 0.1
NS2_start_off                   = TOV

#### system:
SYS_initialize        = TOV+TOV
SYS_initialize_fields = XCTS

#### free data
Free_data_conformal_metric             = flat
Free_data_conformal_Christoffel_symbol = flat
Free_data_conformal_Ricci              = flat
Free_data_trK                          = zero

#### ADM:
ADM_constraints_method     = from_scratch
ADM_B1I_form               = inspiral
ADM_compute_adm_Kuu_method = use_AIJ

#### stress energy tensor:
Tij_NS_decomposition = XCTS
Tij_NS_gConf         = general

#-----------------------------------------------------------------------#
# Settings:                                                             #
#-----------------------------------------------------------------------#

#### checkpoint:
checkpoint_every = 0h

#### basics:
Derivative_Method               = Spectral
Interpolation_Method            = Spectral
Fourier_Transformation_Method   = RFT
dF/du_for_Newton_Method         = Spectral

#-----------------------------------------------------------------------#
# Grid and Geometry:                                                    #
#-----------------------------------------------------------------------#

#### grid:
grid_kind                     = SplitCubedSpherical(NS+NS)
grid_set_NS1                  = left
grid_set_NS2                  = right
grid_NS1_central_box_length   = auto # ~ 20% of the TOV isotropic radius
grid_NS2_central_box_length   = auto # ~ 20% of the TOV isotropic radius
grid_outermost_radius         = 1E5
grid_verbose                  = no

#### resolutions:
n_a = 10(x80)->12(x60)->14(x50)->16(x50)->18(x50)
n_b = 10(x80)->12(x60)->14(x50)->16(x50)->18(x50)
n_c = 12(x80)->14(x60)->16(x50)->18(x50)->20(x50)

grid_SplitCS_max_n_a         = 40
grid_SplitCS_max_n_b         = 40
grid_SplitCS_max_n_c         = 40

#-----------------------------------------------------------------------#
# Equations and Solve:                                                  #
#-----------------------------------------------------------------------#

#### what and where to solve:
Eq_type                         = Elliptic
Eq_elliptic_test                = no

Eq_phi1                         = XCTS_curve_Type3_DDM, NS1
Eq_alias_field_phi1             = phi
Eq_phi2                         = XCTS_curve_Type3_DDM, NS2
Eq_alias_field_phi2             = phi
Eq_psi                          = XCTS_curve_excision_Type1_DDM, .*
Eq_alphaPsi                     = XCTS_curve_excision_Type2_DDM, .*
Eq_B0_U0                        = XCTS_flat_excision_Type1_DDM, .*
Eq_B0_U1                        = XCTS_flat_excision_Type1_DDM, .*
Eq_B0_U2                        = XCTS_flat_excision_Type1_DDM, .*
Eq_update_method                = relaxed_scheme
Eq_update_weight_phi1           = 0.2
Eq_update_weight_phi2           = 0.2
Eq_update_weight_psi            = 0.2
Eq_update_weight_alphaPsi       = 0.2
Eq_update_weight_B0_U0          = 0.2
Eq_update_weight_B0_U1          = 0.2
Eq_update_weight_B0_U2          = 0.2

#### solve settings:
solve_Order                     = phi1,phi2,psi,alphaPsi,B0_U0,B0_U1,B0_U2
solve_Newton_Update_Weight      = 1.
solve_residual                  = 1E-10
solve_residual_factor           = 1E-5
solve_Max_Iteration             = 1
solve_Max_Newton_Step           = 1
solve_Method                    = DDM_Schur_Complement
solve_UMFPACK_refinement_step   = 0
solve_UMFPACK_size              = 1

#-----------------------------------------------------------------------#
# Print:                                                                #
#-----------------------------------------------------------------------#

#### outputs: (regular expression is supported)
txt_output_0d  = ham,mom,eq_residual

#txt_output_1d  = ^phi, ^psi, ^alphaPsi, ^B0, ^beta, eq_residual,
#	         ham, mom, enthalpy, rho0, ^drho0_D.,
#		 d+phi_D.+, d+psi_D.+, d+alphaPsi_D.+, d+B._U.D.+

## lines separated with ',' and the fix values take place in [0,1]x[0,1]
#txt_output_1d_line = (X,0.5,0.5),(0.5,Y,0.5),(0.5,0.5,Z)

#silo_output_3d = ^phi,^psi,^alphaPsi,^B0,^beta,eq_residual,
#                 ham,mom,(dphi_D0,dphi_D1,dphi_D2),
#                 enthalpy,rho0,(drho0_D0,drho0_D1,drho0_D2),
#                 drho0_D0,drho0_D1,drho0_D2,
#	         (B0_U0,B0_U1,B0_U2),(beta_U0,beta_U1,beta_U2)

#-----------------------------------------------------------------------#
# The End                                                               #
#-----------------------------------------------------------------------#

"""


def parse_cli():
  """
    arg parser
    """
  p = argparse.ArgumentParser(
      description="set parameter file for tabulated EoS.")
  p.add_argument(
      "-f",
      type=str,
      required=True,
      help="""/path/to/tab/eos. The table columns must be in the following format:
  "line n(fm^-3) e(g/cm^3) P(dyne/cm^2)"

  Additionally, the specific enthalpy, i.e, (e+p)/rho0, should start from 1.
  """,
  )
  p.add_argument(
      "-m1",
      type=float,
      required=True,
      help="baryonic mass of NS1(geometric unit)",
  )
  p.add_argument(
      "-m2",
      type=float,
      required=True,
      help="baryonic mass of NS2(geometric unit)",
  )
  p.add_argument(
      "-d",
      type=float,
      required=True,
      help="separation between the centers of each NS(geometric unit)",
  )

  args = p.parse_args()
  return args


def eos(args: dict) -> dict:
  """
    smooth eos, if asked, and find values of enthalpy.
    """
  import os
  import numpy as np
  from scipy.interpolate import interp1d

  p_FACTOR = 1.80162095578956e-39
  rho0_FACTOR = 0.002712069678583313
  e_FACTOR = 1.619216164136643e-18
  tab = {}

  fin = args["f"]

  basename = os.path.basename(args["f"]) + "_smooth.txt"
  dirname = os.path.dirname(args["f"])
  fout = os.path.join("./", basename)

  table = np.loadtxt(fin, comments="#")

  num_dens = table[:, 1]
  eng = table[:, 2]
  press = table[:, 3]

  rho0_geo = num_dens * rho0_FACTOR
  press_geo = press * p_FACTOR
  eng_geo = eng * e_FACTOR
  h = (press_geo + eng_geo) / rho0_geo

  h_new = np.linspace(h[0], h[-1], 4 * h.shape[0])
  print("h_new in [a,b] = ", h_new[0], h_new[-1])

  f_e = interp1d(h, eng_geo, kind="linear")
  f_p = interp1d(h, press_geo, kind="linear")

  e_new = f_e(h_new)
  p_new = f_p(h_new)
  r_new = (p_new + e_new) / h_new

  DataOut = np.column_stack((
      np.arange(0, p_new.shape[0], dtype=int),
      r_new / rho0_FACTOR,
      e_new / e_FACTOR,
      p_new / p_FACTOR,
  ))
  np.savetxt(
      fout,
      DataOut,
      header="line n(fm^-3) e(g/cm^3) P(dyne/cm^2)",
      fmt="%d %.16e %.16e %.16e",
  )

  ind_h = np.where(h_new >= 1)
  h = h_new[ind_h]
  tab["file"] = fout
  tab["floor"] = h[0]
  tab["ceiling"] = h_new[-2]

  return tab


def main(args):
  """
    populate place holders
    """

  tab = eos(args)

  pars = {}
  pars["NSNS_separation"] = args["d"]
  pars["NS1_baryonic_mass"] = args["m1"]
  pars["NS2_baryonic_mass"] = args["m2"]
  pars["NS1_EoS_table_path"] = tab["file"]
  pars["NS2_EoS_table_path"] = tab["file"]
  pars["NS1_EoS_enthalpy_floor"] = tab["floor"]
  pars["NS1_EoS_enthalpy_ceiling"] = tab["ceiling"]
  pars["NS2_EoS_enthalpy_floor"] = tab["floor"]
  pars["NS2_EoS_enthalpy_ceiling"] = tab["ceiling"]

  parfile = g_dummy.format(**pars)

  output = f'NS1_m{args["m1"]:0.2f}_O0.0--NS2_m{args["m2"]:0.2f}_O0.0_eostab--d{args["d"]:0.2f}.par'

  with open(f"{output}", "w") as file:
    for k, v in pars.items():
      file.write(f"# {k} = {v}\n")
    file.write(parfile)

  print(f"NOTE:\n{tab['file']} and {output} file \n"
        "must be in the same directory when Elliptica is invoked.")


if __name__ == "__main__":
  args = parse_cli()
  main(args.__dict__)
