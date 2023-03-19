/*
// Alireza Rashti
// January 2023
*/

#include "nsns_main.h"


/* initial data for NS-NS binary system */
int NS_NS_Binary_Initial_Data(void *vp)
{
  /* if this is a generic ID reader call */
  if (strcmp_i(PgetsEZ("IDR_NSNS_export_id"),"generic"))
    nsns_export_id_generic(vp);
  
  /* otherwise construct initial data */
  else
    construct_initial_data(vp);
  
  return EXIT_SUCCESS;
}

/* main algorithm for construction of ID */
static void construct_initial_data(void *vp)
{
  FUNC_TIC
  
  Physics_T *new_phys = 0,*old_phys = 0;
  int Stop = 0;/* if 1 it stops the main loop */
  Uint iter = 0;
  
  /* setting the default parameters */
  set_default_parameters();
  
  /* main iteration loop */
  while(!Stop)
  {
    printf("{ Outermost iteration %u ...\n",iter);
    
    Stop = update_iteration_params(iter,P_,P_"%s_%ux%ux%u");
    
    new_phys = nsns_initialize_new_physics(old_phys);
    
    write_checkpoint(new_phys,Pgets(P_"my_directory"));
    
    nsns_solve_equation(new_phys);
    
    nsns_analyze(new_phys,Pgeti(P_"resolution_iteration"));
    
    free_physics(old_phys);
    
    old_phys = new_phys;
    
    printf("} Outermost iteration %u ==> Done.\n",iter);
    
    iter++;
  }
  
  free_physics(new_phys);
  UNUSED(vp);
  FUNC_TOC
}

/* default parameters for this project */
static void set_default_parameters(void)
{
  /* general parameters: */
  // =================== */
  
  /* which fields to be saved at checkpoint file */
  Pset_default("checkpoint_save","enthalpy,phi,psi,alphaPsi,"
                                 "B0_U0,B0_U1,B0_U2");
  
  /* how to computer derivative */
  Pset_default("Derivative_method","Spectral");
  
  /* how to computer interpolation */
  Pset_default("Interpolation_method","Spectral");
  
  /* how to compute coefficients of expansion:
  // option:
  // =======
  // RFT : regular Fourier transformation. */
  Pset_default("Fourier_transformation_method","RFT");
  
  /* how to computer functional derivatives used in Jacobian */
  Pset_default("dF/du_for_Newton_method","Spectral");
  
  /* optimize ccs reader with 10 splits; if gives error make it smaller */
  Pset_default("matrix_ccs_reader_split","10");
  /* drop matrix entries if less than this. */
  Pset_default("matrix_ccs_drop_below","1.0E-12");
  /* use long UMFPACK */
  Pset_default("solve_UMFPACK_size","1");
  /* set the portion of cores to be used in schur algorithm.
  // the main use of this param is when the memory is scarce and
  // not all the threads should be used simultaneously otherwise 
  // it is killed by the OS. values are in [0,1]. e.g., 0.5 means
  // use 50% of the max number of available threads. */
  Pset_default("solve_ddm_schur_thread_cap","1");
  
  
  /* NS-NS parameters:
  // ================== */
 
  /* how far are NS-NS */
  Pset_default(P_"separation","50.");
  
  /* how fast NS-NS angular velocity */
  Pset_default(P_"angular_velocity","0.");
  
  /* how fast NS-NS infall velocity */
  Pset_default(P_"infall_velocity","0.");
 
  /* how to start off:
  // options:
  // "parameter_file" : it reads parameter file and initialize.
  // "checkpoint_file": it uses checkpoint file and initialize. */
  Pset_default(P_"start_off","parameter_file");
  
  /* total number of iterations that have been executed */
  Pset_default(P_"iteration","0");
  /* total number of iterations executed for this resolution */
  Pset_default(P_"resolution_iteration","0");
  /* stop the main loop if it is 1. */
  Pset_default(P_"STOP","0");
  
  /* system center of mass. */
  Pset_default(P_"x_CM","0."); 
  Pset_default(P_"y_CM","0."); 
  Pset_default(P_"z_CM","0."); 
  
  /* boost velocity for NS-NS */
  Pset_default(P_"boost_Vx","0."); 
  Pset_default(P_"boost_Vy","0."); 
  Pset_default(P_"boost_Vz","0."); 
  
  /* masses */
  Pset_default(P_"ADM_mass","1.");
  Pset_default(P_"Komar_mass","1.");
  
  /* what to print for properties of NSNS, add and separate with comma */
  Pset_default(P_"NSNS_properties",
    "separation,x_CM,y_CM,z_CM,ADM_mass,Komar_mass,mass_ratio,"
    "angular_velocity,infall_velocity,"
    "Px_ADM,Py_ADM,Pz_ADM,"
    "Jx_ADM,Jy_ADM,Jz_ADM,"
    "binding_energy,Virial_error");
    
  /* how to tune P_ADM:
  // options:[none,adjust(?_CM)], cf. */
  Pset_default(P_"P_ADM_control_method","adjust(x_CM,y_CM)");
  Pset_default(P_"P_ADM_control_update_weight","0.");
  Pset_default(P_"P_ADM_control_tolerance","1E-5");
  Pset_default(P_"P_ADM_control_threshold","1E-1");
  
  /* observer method */
  Pset_default(P_"Observe_ADM_P","S+V,constraint");
  Pset_default(P_"Observe_ADM_J","S+V,constraint");
  Pset_default(P_"Observe_ADM_M","S+V,conformal");
  Pset_default(P_"Observe_Komar_M","S_inf,default");
  
  /* equation */
  /* kind of dF/du to compute before calling elliptic solve.
  // if set "none" it won't prepare in advance. 
  // this is used for optimization.
  // NOTE: to optimize further it's better to start from second order. */
  //Pset_default(P_"dF/du_prepare","J_D0D0,J_D0D1,J_D0D2,"
  //                               "J_D1D1,J_D1D2,J_D2D2,"
  //                               "J_D0,J_D1,J_D2");
  /* NOTE: a fast and analytic method implemented thus no further 
  // dF/du preparation is requried. */
  Pset_default(P_"dF/du_prepare","none");
  
  /* NS1 paramters:
  // ============= */
  
  /* NS masses */
  Pset_default("NS1_baryonic_mass_current","1.");
  Pset_default("NS1_baryonic_mass","1.");
  Pset_default("NS1_ADM_mass","1.");
  Pset_default("NS1_Komar_mass","1.");
  Pset_default("NS1_TOV_ADM_mass","1.");
  Pset_default("NS1_TOV_radius","1.");
  Pset_default("NS1_TOV_compactness","0.");
  Pset_default("NS1_mass_shedding_indicator","1.");
  
  /* NS EoS: */
  Pset_default("NS1_EoS_description","NA");
  
  /* [polytropic,piecewise_polytropic] */
  Pset_default("NS1_EoS_type","NA");
  
  /* unit: [geo] */
  Pset_default("NS1_EoS_unit","NA");
  Pset_default("NS1_EoS_K0","NA");
  Pset_default("NS1_EoS_Gamma","NA");
  Pset_default("NS1_EoS_rho0_th","NA");
  
  /* -> central matters */
  Pset_default("NS1_rho0_center","1E-3");
  Pset_default("NS1_pressure_center","1E-3");
  Pset_default("NS1_energy_density_center","1E-3");
  
  /* geometrical center of NS.
  // NOTE: geometrical center can be different from patch->c. */ 
  Pset_default("NS1_center_x","0."); 
  Pset_default("NS1_center_y","0."); 
  Pset_default("NS1_center_z","0."); 
  Pset_default("NS1_x_CM","0."); 
  Pset_default("NS1_y_CM","0."); 
  Pset_default("NS1_z_CM","0."); 

  /* box length at the center of NS */
  Pset_default("grid_NS1_central_box_length","auto");
  
  /* spin vector to adjust spin for NS */
  Pset_default("NS1_Omega_x","0."); 
  Pset_default("NS1_Omega_y","0."); 
  Pset_default("NS1_Omega_z","0."); 

  /* spin */
  Pset_default("NS1_chi_x","0."); 
  Pset_default("NS1_chi_y","0."); 
  Pset_default("NS1_chi_z","0."); 
  Pset_default("NS1_spin_x","0."); 
  Pset_default("NS1_spin_y","0."); 
  Pset_default("NS1_spin_z","0."); 

  /* ADM momentum */
  Pset_default("NS1_Px_ADM","0."); 
  Pset_default("NS1_Py_ADM","0."); 
  Pset_default("NS1_Pz_ADM","0."); 
  
  /* ADM angular momentum */
  Pset_default("NS1_Jx_ADM","0."); 
  Pset_default("NS1_Jy_ADM","0."); 
  Pset_default("NS1_Jz_ADM","0."); 

  /* what to print for properties of NS, add and separate with comma */
  Pset_default(P_"NS1_properties",
   "center_x,center_y,center_z,x_CM,y_CM,z_CM,"
   "max_radius,min_radius,TOV_radius,"
   "ADM_mass,TOV_ADM_mass,Komar_mass,baryonic_mass_current,"
   "baryonic_mass,mass_shedding_indicator,TOV_compactness,"
   "rho0_center,pressure_center,energy_density_center,"
   "enthalpy_L2_residual,Euler_equation_constant,"
   "Omega_x,Omega_y,Omega_z,chi_x,chi_y,chi_z,"
   "spin_x,spin_y,spin_z,"
   "Px_ADM,Py_ADM,Pz_ADM,"
   "Jx_ADM,Jy_ADM,Jz_ADM");
  
  /* the very first NS approximation
  // options:
  // ========
  // o. TOV (see NS physics) */
  Pset_default("NS1_start_off","TOV"); 
  
  /* max l in Ylm expansion */
  Pset_default("NS1_surface_Ylm_max_l","1"); 
  
  /* NS surface */
  Pset_default("NS1_surface_type","perfect_s2"); 
  Pset_default("NS1_did_NS_surface_change?","1");
  
  /* if new surface relative change exceeds this, NS surface gets updated */
  Pset_default("NS1_surface_change_threshold","0.0");
  
  /* max allowed surface fails */
  Pset_default("NS1_surface_max_fail","100");
  
  /* number of surface fails */
  Pset_default("NS1_surface_num_fail","0");
  
  /* observe method pertinet to NS */
  Pset_default("NS1_Observe_ADM_M","V_obj,default");
  Pset_default("NS1_Observe_Komar_M","S_obj,default");
  Pset_default("NS1_Observe_baryonic_M","V_obj,default");
  Pset_default("NS1_Observe_ADM_P","S_obj,default");
  Pset_default("NS1_Observe_ADM_J","S_obj,default");
  Pset_default("NS1_Observe_CM","V_obj,default");
  Pset_default("NS1_Observe_spin","S_obj,JRP");
  
  /* smooth and polish phi equation close to the surface */
  Pset_default("NS1_Eq_phi_polish","0.1");
  
  /* tune and adjust: */
  /* for option cf star_main */
  Pset_default("NS1_force_balance_equation","adjust(d/dy:Omega)");
  Pset_default("NS1_force_balance_update_weight","0.");
  Pset_default("NS1_adjust_center_method","Taylor_expansion");
  Pset_default("NS1_adjust_center_update_weight","0.");
  Pset_default("NS1_enthalpy_allowed_residual","1E-8");
  Pset_default("NS1_enthalpy_L2_residual","0.");
  
  /* extrapolation of matter fields outside NS:
  // options = [inverse_r2,exp2,poly2,inverse_r2_expmr,
                inverse_r2_expmAr,enthalpy_expmr_phi_inverse_r2]. */
  Pset_default("NS1_extrapolate_matter_fields","inverse_r2_expmAr");
  
  /* Euler eq. constant */
  Pset_default("NS1_Euler_equation_constant","-0.8");
  Pset_default("NS1_Euler_const_update_weight","1.");
  
  /* NS enhtalpy update weight */
  Pset_default("NS1_enthalpy_update_weight","0.");
  
  /* set enthalpy != 1 on surface to 1
  // options: [yes/no]. */
  Pset_default("NS1_enthalpy_neat","yes");
  
  /* root finder pertinent to NS */
  Pset_default("NS1_RootFinder_method","Steepest_Descent");
  Pset_default("NS1_RootFinder_Tolerance","1E-9");
  Pset_default("NS1_RootFinder_Iteration","1000");
  Pset_default("NS1_RootFinder_verbose","no");
  
  /* NS2 paramters:
  // ============= */
  
  /* NS masses */
  Pset_default("NS2_baryonic_mass_current","1.");
  Pset_default("NS2_baryonic_mass","1.");
  Pset_default("NS2_ADM_mass","1.");
  Pset_default("NS2_Komar_mass","1.");
  Pset_default("NS2_TOV_ADM_mass","1.");
  Pset_default("NS2_TOV_radius","1.");
  Pset_default("NS2_TOV_compactness","0.");
  Pset_default("NS2_mass_shedding_indicator","1.");
  
  /* NS EoS: */
  Pset_default("NS2_EoS_description","NA");
  
  /* [polytropic,piecewise_polytropic] */
  Pset_default("NS2_EoS_type","NA");
  
  /* unit: [geo] */
  Pset_default("NS2_EoS_unit","NA");
  Pset_default("NS2_EoS_K0","NA");
  Pset_default("NS2_EoS_Gamma","NA");
  Pset_default("NS2_EoS_rho0_th","NA");
  
  /* -> central matters */
  Pset_default("NS2_rho0_center","1E-3");
  Pset_default("NS2_pressure_center","1E-3");
  Pset_default("NS2_energy_density_center","1E-3");
  
  /* geometrical center of NS.
  // NOTE: geometrical center can be different from patch->c. */ 
  Pset_default("NS2_center_x","0."); 
  Pset_default("NS2_center_y","0."); 
  Pset_default("NS2_center_z","0."); 
  Pset_default("NS2_x_CM","0."); 
  Pset_default("NS2_y_CM","0."); 
  Pset_default("NS2_z_CM","0."); 
  
  /* box length at the center of NS */
  Pset_default("grid_NS2_central_box_length","auto");
  
  /* spin vector to adjust spin for NS */
  Pset_default("NS2_Omega_x","0."); 
  Pset_default("NS2_Omega_y","0."); 
  Pset_default("NS2_Omega_z","0."); 
  
  /* spin */
  Pset_default("NS2_chi_x","0."); 
  Pset_default("NS2_chi_y","0."); 
  Pset_default("NS2_chi_z","0."); 
  Pset_default("NS2_spin_x","0."); 
  Pset_default("NS2_spin_y","0."); 
  Pset_default("NS2_spin_z","0."); 
  
  /* ADM momentum */
  Pset_default("NS2_Px_ADM","0."); 
  Pset_default("NS2_Py_ADM","0."); 
  Pset_default("NS2_Pz_ADM","0."); 
  
  /* ADM angular momentum */
  Pset_default("NS2_Jx_ADM","0."); 
  Pset_default("NS2_Jy_ADM","0."); 
  Pset_default("NS2_Jz_ADM","0."); 
  
  /* what to print for properties of NS, add and separate with comma */
  Pset_default(P_"NS2_properties",
   "center_x,center_y,center_z,x_CM,y_CM,z_CM,"
   "max_radius,min_radius,TOV_radius,"
   "ADM_mass,TOV_ADM_mass,Komar_mass,baryonic_mass_current,"
   "baryonic_mass,mass_shedding_indicator,TOV_compactness,"
   "rho0_center,pressure_center,energy_density_center,"
   "enthalpy_L2_residual,Euler_equation_constant,"
   "Omega_x,Omega_y,Omega_z,chi_x,chi_y,chi_z,"
   "spin_x,spin_y,spin_z,"
   "Px_ADM,Py_ADM,Pz_ADM,"
   "Jx_ADM,Jy_ADM,Jz_ADM");
  
  /* the very first NS approximation
  // options:
  // ========
  // o. TOV (see NS physics) */
  Pset_default("NS2_start_off","TOV"); 
  
  /* max l in Ylm expansion */
  Pset_default("NS2_surface_Ylm_max_l","1"); 
  
  /* NS surface */
  Pset_default("NS2_surface_type","perfect_s2"); 
  Pset_default("NS2_did_NS_surface_change?","1");
  
  /* if new surface relative change exceeds this, NS surface gets updated */
  Pset_default("NS2_surface_change_threshold","0.0");
  
  /* max allowed surface fails */
  Pset_default("NS2_surface_max_fail","100");
  
  /* number of surface fails */
  Pset_default("NS2_surface_num_fail","0");
  
  /* observe method pertinet to NS */
  Pset_default("NS2_Observe_ADM_M","V_obj,default");
  Pset_default("NS2_Observe_Komar_M","S_obj,default");
  Pset_default("NS2_Observe_baryonic_M","V_obj,default");
  Pset_default("NS2_Observe_ADM_P","S_obj,default");
  Pset_default("NS2_Observe_ADM_J","S_obj,default");
  Pset_default("NS2_Observe_CM","V_obj,default");
  Pset_default("NS2_Observe_spin","S_obj,JRP");
  
  /* smooth and polish phi equation close to the surface */
  Pset_default("NS2_Eq_phi_polish","0.1");
  
  /* tune and adjust: */
  /* for option cf star_main */
  Pset_default("NS2_force_balance_equation","adjust(d/dy:Omega)");
  Pset_default("NS2_force_balance_update_weight","0.");
  Pset_default("NS2_adjust_center_method","Taylor_expansion");
  Pset_default("NS2_adjust_center_update_weight","0.");
  Pset_default("NS2_enthalpy_allowed_residual","1E-8");
  Pset_default("NS2_enthalpy_L2_residual","0.");
  
  /* extrapolation of matter fields outside NS:
  // options = [inverse_r2,exp2,poly2,inverse_r2_expmr,
                inverse_r2_expmAr,enthalpy_expmr_phi_inverse_r2]. */
  Pset_default("NS2_extrapolate_matter_fields","inverse_r2_expmAr");
  
  /* Euler eq. constant */
  Pset_default("NS2_Euler_equation_constant","-0.8");
  Pset_default("NS2_Euler_const_update_weight","1.");
  
  /* NS enhtalpy update weight */
  Pset_default("NS2_enthalpy_update_weight","0.");
  
  /* set enthalpy != 1 on surface to 1
  // options: [yes/no]. */
  Pset_default("NS2_enthalpy_neat","yes");
  
  /* root finder pertinent to NS */
  Pset_default("NS2_RootFinder_method","Steepest_Descent");
  Pset_default("NS2_RootFinder_Tolerance","1E-9");
  Pset_default("NS2_RootFinder_Iteration","1000");
  Pset_default("NS2_RootFinder_verbose","no");
  
  /* use Kepler's law to set angular velocity */
  if (Pcmps(P_"angular_velocity","auto"))
  {
    const double m1 = Pgetd("NS1_baryonic_mass");
    const double m2 = Pgetd("NS2_baryonic_mass");
    
    const double r  = Pgetd(P_"separation");
    const double O  = sqrt((m1+m2)/pow(r,3.));
    
    Psetd(P_"angular_velocity",O);
  }
  
  /* set NS1 and NS2 centers */
  if (strstr(Pgets("grid_set_NS1"),"right"))
  {
    const double S = Pgetd(P_"separation");
    
    /* NS1 center in +y */
    Psetd("NS1_center_x",0.);
    Psetd("NS1_center_y",S/2.);
    Psetd("NS1_center_z",0.);
    
    /* NS2 center in -y */
    Psetd("NS2_center_x",0.);
    Psetd("NS2_center_y",-S/2.);
    Psetd("NS2_center_z",0.);
  }
  else
  {
    const double S = Pgetd(P_"separation");

    /* NS1 center in -y */
    Psetd("NS1_center_x",0.);
    Psetd("NS1_center_y",-S/2.);
    Psetd("NS1_center_z",0.);
    
    /* NS2 center in +y */
    Psetd("NS2_center_x",0.);
    Psetd("NS2_center_y",S/2.);
    Psetd("NS2_center_z",0.);
    
  }
  
  /* set center of mass */
  {
    const double bar_mass1 = Pgetd("NS1_baryonic_mass");
    const double bar_mass2 = Pgetd("NS2_baryonic_mass");
    const double y_CM = (bar_mass1*Pgetd("NS1_center_y") +
                         bar_mass2*Pgetd("NS2_center_y"))/(bar_mass1+bar_mass2);
   
    Psetd(P_"x_CM",0.);
    Psetd(P_"y_CM",y_CM);
    Psetd(P_"z_CM",0.);
  }
  
}
