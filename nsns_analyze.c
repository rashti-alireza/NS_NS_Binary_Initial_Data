/*
// Alireza Rashti
// January 2023
*/

/* analyzing initial data such as mass, momentum, constraints etc.  */

#include "nsns_analyze.h"

/* analyzing physics properties, constraints etc */
void nsns_analyze(Physics_T *const phys,const int iteration)
{
  if (!phys) return;

  FUNC_TIC
  
  const char *const properties_file_name = P_"properties.txt";
  FILE *file = 0;
  char str[MAX_STR_LEN];
   
  /* compute properties and constraints */ 
  physics(phys,ADM_COMPUTE_CONSTRAINTS);
  
  /* compute various properties */
  compute_properties(phys);
  
  /* open properties file in "my_directory" and save */
  sprintf(str,"%s/%s",Pgets(P_"my_directory"),properties_file_name);
  file = Fopen(str,"w");
  nsns_print_physical_system_properties(phys,file,iteration,0);
  Fclose(file);

  /* open properties file in "Diagnostics" and save */
  sprintf(str,"%s/%s",Pgets(P_"Diagnostics"),properties_file_name);
  file = Fopen(str,"a");
  nsns_print_physical_system_properties(phys,file,iteration,0);
  Fclose(file);
  
  /* prints */
  print_fields_0D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_1D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_2D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_3D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  
  FUNC_TOC
}

/* print physical system properties such as mass, spin etc in the given
// file, if pr_screen is 1, it also prints in stdout */
void nsns_print_physical_system_properties(Physics_T *const phys,
                                          FILE *const file,
                                          const int iteration,
                                          const int pr_screen)
{
  Physics_T *const ns1 = init_physics(phys,NS1);
  Physics_T *const ns2 = init_physics(phys,NS2);

  if (pr_screen)
  {
    printf(Pretty0"iteration = %d:\n",iteration);
  }
  fprintf(file,"%s\n",LINE_STR);
  fprintf(file,"# iteration = %d\n",iteration);
  fprintf(file,"\n");
  
  star_print_properties(ns1,Pgets(P_"NS1_properties"),file,pr_screen);
  star_print_properties(ns2,Pgets(P_"NS2_properties"),file,pr_screen);
  sys_print_properties(phys,Pgets(P_"NSNS_properties"),file,pr_screen);
  
  free_physics(ns1);
  free_physics(ns2);
}

/* compute variety of properties.
// NOTE: order of parameter calculations matter. 
// NOTE: if there is a confusion between target params and current
//       params, "current" suffix added to the latter. */
static void compute_properties(Physics_T *const phys/* nsns */)
{
  Physics_T *const ns1 = init_physics(phys,NS1);
  Physics_T *const ns2 = init_physics(phys,NS2);
  TOV_T *tov        = 0;
  
  const double x_CM = Pgetd(P_"x_CM");
  const double y_CM = Pgetd(P_"y_CM");
  const double z_CM = Pgetd(P_"z_CM");
  double p[3]  = {0.};
  double j[3]  = {0.};
  double s[3]  = {0.};
  double cm[3] = {0.};
  double m     = 0.;
  
  /* NS1: */
  observe(ns1,"ADM(M)",Pgets("NS1_Observe_ADM_M"),&m);
  Psetd("NS1_ADM_mass",m);
  
  observe(ns1,"Komar(M)",Pgets("NS1_Observe_Komar_M"),&m);
  Psetd("NS1_Komar_mass",m);
  
  observe(ns1,"Baryonic(M)",Pgets("NS1_Observe_baryonic_M"),&m);
  Psetd("NS1_baryonic_mass_current",m);
  
  tov = TOV_init();
  tov->exit_if_error = 0;
  tov->phys  = ns1;
  tov->bar_m = Pgetd("NS1_baryonic_mass_current");
  tov = TOV_solution(tov);
  if (tov->status == 0)
  {
    Psetd("NS1_TOV_ADM_mass",tov->ADM_m);
    /* Note: compactness = adm_mass/Schwarzschild_radius 
      (not isotropic radius) */
    Psetd("NS1_TOV_compactness",tov->ADM_m/tov->r[tov->N-1]);
    Psetd("NS1_TOV_radius",tov->rbar[tov->N-1]);
  }
  TOV_free(tov);
  
  Psetd("NS1_mass_shedding_indicator",star_NS_mass_shedding_indicator(ns1));
  
  observe(ns1,"CM",Pgets("NS1_Observe_CM"),cm);
  Psetd("NS1_x_CM",cm[0]+x_CM);
  Psetd("NS1_y_CM",cm[1]+y_CM);
  Psetd("NS1_z_CM",cm[2]+z_CM);

  observe(ns1,"ADM(P)",Pgets("NS1_Observe_ADM_P"),p);
  Psetd("NS1_Px_ADM",p[0]);
  Psetd("NS1_Py_ADM",p[1]);
  Psetd("NS1_Pz_ADM",p[2]);

  observe(ns1,"ADM(J)",Pgets("NS1_Observe_ADM_J"),j);
  Psetd("NS1_Jx_ADM",j[0]);
  Psetd("NS1_Jy_ADM",j[1]);
  Psetd("NS1_Jz_ADM",j[2]);
  
  observe(ns1,"spin",Pgets("NS1_Observe_spin"),s);
  Psetd("NS1_Spin_x",s[0]);
  Psetd("NS1_Spin_y",s[1]);
  Psetd("NS1_Spin_z",s[2]);
  
  m = Pgetd("NS1_adm_mass");
  Psetd("NS1_chi_x",s[0]/Pow2(m));
  Psetd("NS1_chi_y",s[1]/Pow2(m));
  Psetd("NS1_chi_z",s[2]/Pow2(m));
  
  /* NS2: */
  observe(ns2,"ADM(M)",Pgets("NS2_Observe_ADM_M"),&m);
  Psetd("NS2_ADM_mass",m);
  
  observe(ns2,"Komar(M)",Pgets("NS2_Observe_Komar_M"),&m);
  Psetd("NS2_Komar_mass",m);
  
  observe(ns2,"Baryonic(M)",Pgets("NS2_Observe_baryonic_M"),&m);
  Psetd("NS2_baryonic_mass_current",m);
  
  tov = TOV_init();
  tov->exit_if_error = 0;
  tov->phys  = ns2;
  tov->bar_m = Pgetd("NS2_baryonic_mass_current");
  tov = TOV_solution(tov);
  if (tov->status == 0)
  {
    Psetd("NS2_TOV_ADM_mass",tov->ADM_m);
    /* Note: compactness = adm_mass/Schwarzschild_radius 
      (not isotropic radius) */
    Psetd("NS2_TOV_compactness",tov->ADM_m/tov->r[tov->N-1]);
    Psetd("NS2_TOV_radius",tov->rbar[tov->N-1]);
  }
  TOV_free(tov);
  
  Psetd("NS2_mass_shedding_indicator",star_NS_mass_shedding_indicator(ns2));
  
  observe(ns2,"CM",Pgets("NS2_Observe_CM"),cm);
  Psetd("NS2_x_CM",cm[0]+x_CM);
  Psetd("NS2_y_CM",cm[1]+y_CM);
  Psetd("NS2_z_CM",cm[2]+z_CM);

  observe(ns2,"ADM(P)",Pgets("NS2_Observe_ADM_P"),p);
  Psetd("NS2_Px_ADM",p[0]);
  Psetd("NS2_Py_ADM",p[1]);
  Psetd("NS2_Pz_ADM",p[2]);

  observe(ns2,"ADM(J)",Pgets("NS2_Observe_ADM_J"),j);
  Psetd("NS2_Jx_ADM",j[0]);
  Psetd("NS2_Jy_ADM",j[1]);
  Psetd("NS2_Jz_ADM",j[2]);
  
  observe(ns2,"spin",Pgets("NS2_Observe_spin"),s);
  Psetd("NS2_Spin_x",s[0]);
  Psetd("NS2_Spin_y",s[1]);
  Psetd("NS2_Spin_z",s[2]);
  
  m = Pgetd("NS2_adm_mass");
  Psetd("NS2_chi_x",s[0]/Pow2(m));
  Psetd("NS2_chi_y",s[1]/Pow2(m));
  Psetd("NS2_chi_z",s[2]/Pow2(m));
  
  
  /* NSNS: */
  observe(phys,"ADM(M)",Pgets(P_"Observe_ADM_M"),&m);
  Psetd(P_"adm_mass",m);
  
  observe(phys,"Komar(M)",Pgets(P_"Observe_Komar_M"),&m);
  Psetd(P_"Komar_mass",m);
  
  observe(phys,"ADM(P)",Pgets(P_"Observe_ADM_P"),p);
  Psetd(P_"Px_ADM",p[0]);
  Psetd(P_"Py_ADM",p[1]);
  Psetd(P_"Pz_ADM",p[2]);

  observe(phys,"ADM(J)",Pgets(P_"Observe_ADM_J"),j);
  Psetd(P_"Jx_ADM",j[0]);
  Psetd(P_"Jy_ADM",j[1]);
  Psetd(P_"Jz_ADM",j[2]);

  /* mass ratio (q >= 1.) */
  double q = Pgetd("NS2_TOV_ADM_mass")/Pgetd("NS1_TOV_ADM_mass");
  Psetd(P_"mass_ratio",q > 1. ? q : 1./q );

  /* binding energy */
  double bin_e = Pgetd(P_"adm_mass") -
    (Pgetd("NS2_TOV_ADM_mass") + Pgetd("NS1_TOV_ADM_mass"));
  Psetd(P_"binding_energy",bin_e);

  /* virial error */
  double v_e = fabs(1.-Pgetd(P_"adm_mass")/Pgetd(P_"komar_mass"));
  Psetd(P_"virial_error",v_e);
  
  /* number of orbits (lowest order PN) */
  if (0)/* too much inaccurate */
  {
    double m1    = Pgetd("NS1_adm_mass");
    double m2    = Pgetd("NS2_adm_mass");
    double nu    = m1*m2/Pow2(m1+m2);
    double omega = Pgetd(P_"angular_velocity");
    double m_tot = Pgetd(P_"ADM_mass");
    double N_orb = pow(m_tot*omega,-5./3.)/(32.*nu)/(2.*M_PI);
    /* initially some vars might be off or negative like total adm_mass */
    N_orb = (isfinite(N_orb) && N_orb > 0. ? N_orb : 0.);
    Psetd(P_"number_of_orbits_1PN",N_orb);
  }
  
  free_physics(ns1);
  free_physics(ns2);
}
