/*
 * kmaplib_models.cpp
 * 
 * This file contains all the kinetic modeling functions used to calculate time 
 * activity curves (TAC) and sensitivity functions for different kinetic models.
 */

#include "kmaplib.h"
#include "mex.h"

//------------------------------------------------------------------------------
// tac_eval
//------------------------------------------------------------------------------
// Evaluates the time activity curve (TAC) for the provided kinetic model parameters.
void tac_eval(double *p, void *param, double *tac)
{ 
   KMODEL_T *par;
   par = (KMODEL_T *) param;
   (*(par->tacfunc))(p, par->dk, par->scant, par->td, par->cp, par->wb, 
                     par->num_frm, par->num_vox, tac);
}

//------------------------------------------------------------------------------
// jac_eval
//------------------------------------------------------------------------------
// Evaluates the Jacobian (sensitivity matrix) for the provided kinetic model parameters.
void jac_eval(double *p, void *param, double *tac, int *psens, double *jac)
{
   KMODEL_T *par;
   par = (KMODEL_T *) param;
   (*(par->jacfunc))(p, par->dk, par->scant, par->td, par->cp, par->wb, 
                     par->num_frm, par->num_vox, tac, psens, jac);
}

//------------------------------------------------------------------------------
// kconv_2tcm_tac
//------------------------------------------------------------------------------
// Calculates the time activity curve (TAC) using the two-tissue compartment model.
// This model uses six parameters: vb, k1, k2, k3, k4, and t_delay.
void kconv_2tcm_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct)
{
   int      i, j;
   double   vb, k1, k2, k3, k4, t_delay;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double  *c_a1, *c_a2, *c_f, *c_b, *c_t;
   double  *cp_delay, *wb_delay;  
   int      num_time;

   // Memory allocation for intermediate calculations
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a1 = (double*) malloc(sizeof(double)*num_time);
   c_a2 = (double*) malloc(sizeof(double)*num_time);
   c_f  = (double*) malloc(sizeof(double)*num_time);
   c_b  = (double*) malloc(sizeof(double)*num_time);
   c_t  = (double*) malloc(sizeof(double)*num_time);
   cp_delay = (double*) malloc(sizeof(double)*num_time);
   wb_delay = (double*) malloc(sizeof(double)*num_time);

   for (j=0; j<num_vox; j++) {
      // Parameter transformations
      vb = p[0+j*6];
      k1 = p[1+j*6];
      k2 = p[2+j*6];
      k3 = p[3+j*6];
      k4 = p[4+j*6];
      t_delay = p[5+j*6];
      k234 = k2+k3+k4;
      d = sqrt(k234*k234 - 4*k2*k4);
      a1 = (k234 - d) / 2;
      a2 = (k234 + d) / 2;

      // get delayed cp and wb
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);
      time_delay_tac(cp, num_time, t_delay, td, cp_delay);

      // Convolution with exponential functions
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp_delay, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp_delay, num_time, td, c_a2);

      // Compute compartmental TACs
      if (d==0) d = 1e9;
      f1 = k1/d * (k4 - a1);
      f2 = k1/d * (a2 - k4);
      b1 = k1/d * k3;
      b2 = -b1;

      for (i=0; i<num_time; i++) {
         c_f[i] = f1*c_a1[i] + f2*c_a2[i];
         c_b[i] = b1*c_a1[i] + b2*c_a2[i];
         c_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * wb_delay[i];
      }

      // Frame-averaged activity
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
   }

   // Free allocated memory
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);
   free(c_t);
   free(cp_delay);
   free(wb_delay);
}

//------------------------------------------------------------------------------
// kconv_2tcm_jac
//------------------------------------------------------------------------------
// Calculate the time activity curves and the sensitivity functions 
// (Jacobian) for the two-tissue compartment model.
void kconv_2tcm_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
{
   int      i, j;
   double   vb, k1, k2, k3, k4, t_delay;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double  *c_a1, *c_a2, *c_f, *c_b, *c_t, *s_t;
   double  *wb_delay, *cp_delay;       
   double  *cp_grad_delay, *wb_grad_delay;                          
   int      num_time;
   int      num_par;
  
   // Determine the number of sensitive parameters
   num_par = 0;
   for (i = 0; i < 6; i++) {
      if (psens[i] == 1)
         ++num_par;
   }
  
   // Allocate memory for intermediate calculations
   num_time = (int) (scant[2*num_frm-1] / td);
   c_a1 = (double*) malloc(sizeof(double) * num_time);
   c_a2 = (double*) malloc(sizeof(double) * num_time);
   c_f  = (double*) malloc(sizeof(double) * num_time);
   c_b  = (double*) malloc(sizeof(double) * num_time);

   // shift cp, wb for time delay
   wb_delay  = (double*) malloc(sizeof(double)*num_time);
   cp_delay  = (double*) malloc(sizeof(double)*num_time);
	
   c_t  = (double*) malloc(sizeof(double) * num_time);
   s_t  = (double*) malloc(sizeof(double) * num_time * num_par);
   cp_grad_delay  = (double*) malloc(sizeof(double)*num_time);
   wb_grad_delay  = (double*) malloc(sizeof(double)*num_time);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {
	
      // Transform kinetic parameters into exponential parameters
      vb = p[0 + j * 6];
      k1 = p[1 + j * 6];
      k2 = p[2 + j * 6];
      k3 = p[3 + j * 6];
      k4 = p[4 + j * 6];
      t_delay = p[5+j*6];

      k234 = k2 + k3 + k4;
      d = sqrt(k234 * k234 - 4 * k2 * k4);
      a1 = (k234 - d) / 2;
      a2 = (k234 + d) / 2;
      
      // cp, wb with time delay
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);
      time_delay_tac(cp, num_time, t_delay, td, cp_delay);

      // exp components without time delay
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp_delay, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp_delay, num_time, td, c_a2);

      // Calculate concentrations
      if (d == 0) d = 1e9;
      f1 = k1 / d * (k4 - a1);
      f2 = k1 / d * (a2 - k4);
      b1 = k1 / d * k3;
      b2 = -b1;	
      for (i = 0; i < num_time; i++) {
         c_f[i] = f1 * c_a1[i] + f2 * c_a2[i];
         c_b[i] = b1 * c_a1[i] + b2 * c_a2[i];
         c_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * wb_delay[i];
      }

      // Sensitivity functions (Jacobian) computation
      if (psens[0] == 1) { // Sensitivity wrt vb
         for (i = 0; i < num_time; i++)
            s_t[i] = - (c_f[i] + c_b[i]) + wb_delay[i];
         s_t += num_time;
      }
      if (psens[1] == 1 || psens[2] == 1) { 
         f1 = 1 / d * (k4 + k3 - a1);
         f2 = 1 / d * (a2 - k4 - k3);
      }
      if (psens[1] == 1) { // Sensitivity wrt k1
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[2] == 1 || psens[3] == 1) { 
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a2);
      }
      if (psens[2] == 1) { // Sensitivity wrt k2
         for (i = 0; i < num_time; i++)
            s_t[i] = -(1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[3] == 1 || psens[4] == 1) { 
         f1 = 1 / d * (a1 + a2 - k3 - k4);
         f2 = -f1;
      }
      if (psens[3] == 1) { // Sensitivity wrt k3
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[4] == 1) { // Sensitivity wrt k4
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a2);
         for (i = 0; i < num_time; i++)
            s_t[i] = -(1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[5] == 1) { // wrt time delay
         time_delay_jac(cp_delay, num_time, td, cp_grad_delay);
         time_delay_jac(wb_delay, num_time, td, wb_grad_delay);
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, cp_grad_delay, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, cp_grad_delay, num_time, td, c_a2);
         f1 = k1 / d * (k4 - a1);
         f2 = k1 / d * (a2 - k4);
         b1 = k1 / d * k3;
         b2 = -b1;	
	 for (i=0; i<num_time; i++){ 
            c_f[i] = f1 * c_a1[i] + f2 * c_a2[i];
            c_b[i] = b1 * c_a1[i] + b2 * c_a2[i];
            s_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb *wb_grad_delay[i];
            s_t[i] = - s_t[i];
         } 
	 s_t += num_time;
      }

      s_t -= num_time * num_par;

      // Frame-averaged activity and sensitivity
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
      frame(scant, td, s_t, num_frm, num_par, st + j * num_frm * num_par);
   }
  	
   // Free allocated memory
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);
   free(cp_delay);
   free(wb_delay);
   free(cp_grad_delay);
   free(wb_grad_delay);
   free(c_t);
   free(s_t);
}

//------------------------------------------------------------------------------
// kconv_1tcm_tac
//------------------------------------------------------------------------------
// Calculates the time activity curve (TAC) using the one-tissue compartment model.
// This model uses four parameters: vb, k1, k2, and t_delay.
void kconv_1tcm_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct)
{
   int      i, j;
   double   vb, k1, k2, t_delay;
   double   a;
   double  *c_a, *c_t;
   double  *cp_delay, *wb_delay;      
   int      num_time;

   // Memory allocation
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a = (double*) malloc(sizeof(double)*num_time);
   c_t = (double*) malloc(sizeof(double)*num_time);

   cp_delay = (double*) malloc(sizeof(double)*num_time);
   wb_delay = (double*) malloc(sizeof(double)*num_time);

   for (j = 0; j < num_vox; j++) {
      // Assign parameters
      vb = p[0 + j * 4];
      k1 = p[1 + j * 4];
      k2 = p[2 + j * 4];
      t_delay = p[3+j*4];

      time_delay_tac(cp, num_time, t_delay, td, cp_delay);
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);
      // Convolution with exponential functions
      a = k2 + dk;
      kconv_exp(k1, a, cp_delay, num_time, td, c_a);

      // Compute tissue concentration
      for (i = 0; i < num_time; i++)
         c_t[i] = (1-vb)*c_a[i] + vb*wb_delay[i];

      // Frame-averaged activity
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
   }

   // Free allocated memory
   free(c_a);
   free(c_t);
   free(cp_delay);
   free(wb_delay);
}

//------------------------------------------------------------------------------
// kconv_1tcm_jac
//------------------------------------------------------------------------------
// Calculates the time activity curves and the sensitivity functions for the 
// one-tissue kinetic model.
void kconv_1tcm_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
{
   int      i, j;
   double   vb, k1, k2, t_delay;
   double   a;
   double  *c_a, *c_t, *s_t;
   double  *cp_delay, *wb_delay;  
   double  *cp_grad_delay, *wb_grad_delay;                         
   int      num_time;
   int      num_par;

   // Determine the number of sensitive parameters
   num_par = 0;
   for (i = 0; i < 4; i++) {
      if (psens[i] == 1) ++num_par;
   }

   // Allocate memory for intermediate calculations
   num_time = (int) (scant[2*num_frm-1] / td);
   c_a = (double*) malloc(sizeof(double) * num_time);
   c_t = (double*) malloc(sizeof(double) * num_time);
   s_t = (double*) malloc(sizeof(double) * num_time * num_par);

   cp_delay = (double*) malloc(sizeof(double)*num_time);
   wb_delay = (double*) malloc(sizeof(double)*num_time);
   cp_grad_delay  = (double*) malloc(sizeof(double)*num_time);
   wb_grad_delay  = (double*) malloc(sizeof(double)*num_time);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {

      // Assign parameters
      vb = p[0 + j * 4];
      k1 = p[1 + j * 4];
      k2 = p[2 + j * 4];
      t_delay = p[3 + j * 4];

      // shift cp, wb for time delay
      time_delay_tac(cp, num_time, t_delay, td, cp_delay);
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);

      // Calculate the exponential component
      a = k2 + dk;
      kconv_exp(k1, a, cp_delay, num_time, td, c_a);

      // Calculate the time activity curve
      for (i = 0; i < num_time; i++){
         c_t[i] = (1 - vb) * c_a[i] + vb * wb_delay[i];
      }

      // Calculate sensitivity functions
      if (psens[0] == 1) { // Sensitivity wrt vb
         for (i = 0; i < num_time; i++)
            s_t[i] = -c_a[i] + wb_delay[i];
         s_t += num_time;
      }
      if (psens[1] == 1) { // Sensitivity wrt k1
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * c_a[i] / k1;
         s_t += num_time;
      }
      if (psens[2] == 1) { // Sensitivity wrt k2
         kconv_exp(1, a, c_a, num_time, td, s_t);
         for (i = 0; i < num_time; i++)
            s_t[i] *= -(1 - vb);
         s_t += num_time;
      }
      if (psens[3] == 1) { // wrt time delay
         time_delay_jac(cp_delay, num_time, td, cp_grad_delay);
         time_delay_jac(wb_delay, num_time, td, wb_grad_delay);
         a = k2 + dk;
         kconv_exp(k1, a, cp_grad_delay, num_time, td, c_a);
	      for (i=0; i<num_time; i++){ 
            s_t[i] = (1 - vb) * c_a[i] + vb * wb_grad_delay[i];
            s_t[i] = - s_t[i];
         } 
	 s_t += num_time;
      }
      s_t -= num_time * num_par;

      // Average the sensitivity functions over the frame duration
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
      frame(scant, td, s_t, num_frm, num_par, st + j * num_frm * num_par);
   }

   // Free allocated memory
   free(c_a);
   free(c_t);
   free(s_t);
   free(cp_delay);
   free(wb_delay);
   free(cp_grad_delay);
   free(wb_grad_delay);
}

//------------------------------------------------------------------------------
// kconv_liver_tac
//------------------------------------------------------------------------------
// Calculates the time activity curve using the two-tissue kinetic model with 
// a dual-blood input function model for the liver. 
/* More details of the model are referred to the references:
   [1] Wang GB, Corwin MT, Olson KA, Badawi RD, Sarkar S. Dynamic PET of human 
   liver inflammation: impact of kinetic modeling with optimization-derived dual- 
    blood input function. Physics in Medicine and Biology, 63(15): 155004 (14pp),
    2018.
   [2] Zuo Y, Sarkar S, Corwin MT, Olson K, Badawi RD, Wang GB. Structural and 
   practical identifiability of dual-input kinetic modeling in dynamic PET of liver 
   inflammation. Physics in Medicine and Biology,  64(17): 175023 (18pp),  2019.
*/
void kconv_liver_tac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct)
{
   int      i, j;
   double   vb, k1, k2, k3, k4, ka, fa, t_delay;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double   *c_pv, *cp, *c_a1, *c_a2, *c_f, *c_b, *c_t;
   double   *ca_delay, *wb_delay;              
   int      num_time;

   // Allocate memory for intermediate calculations
   num_time = (int) (scant[2*num_frm-1] / td);
   c_pv = (double*) malloc(sizeof(double) * num_time);
   cp   = (double*) malloc(sizeof(double) * num_time);
   c_a1 = (double*) malloc(sizeof(double) * num_time);
   c_a2 = (double*) malloc(sizeof(double) * num_time);
   c_f  = (double*) malloc(sizeof(double) * num_time);
   c_b  = (double*) malloc(sizeof(double) * num_time);
   c_t  = (double*) malloc(sizeof(double) * num_time);
   ca_delay = (double*) malloc(sizeof(double)*num_time);
   wb_delay = (double*) malloc(sizeof(double)*num_time);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {

      // Transform kinetic parameters into exponential parameters
      vb = p[0+j*8];
      k1 = p[1+j*8];
      k2 = p[2+j*8];
      k3 = p[3+j*8];
      k4 = p[4+j*8];
      ka = p[5+j*8];
      fa = p[6+j*8];
      t_delay = p[7 + j*8];

      k234 = k2 + k3 + k4;
      d = sqrt(k234 * k234 - 4 * k2 * k4);
      a1 = (k234 - d) / 2;
      a2 = (k234 + d) / 2;
      // get delayed cp/wb
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);
      time_delay_tac(ca, num_time, t_delay, td, ca_delay);

      // Calculate dispersion and input concentration
      tmp = ka + dk;
      kconv_exp(ka, tmp, ca_delay, num_time, td, c_pv);
      for (i = 0; i < num_time; i++) {
         cp[i] = (1 - fa) * c_pv[i] + fa * ca_delay[i];
      }

      // Calculate the exponential components
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2);

      // Calculate the concentrations
      if (d == 0) d = 1e9;
      f1 = k1 / d * (k4 - a1);
      f2 = k1 / d * (a2 - k4);
      b1 = k1 / d * k3;
      b2 = -b1;
      for (i = 0; i < num_time; i++) {
         c_f[i] = f1 * c_a1[i] + f2 * c_a2[i];
         c_b[i] = b1 * c_a1[i] + b2 * c_a2[i];
         c_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * cp[i];
      }

      // Average the values over frame duration
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
   }

   // Free allocated memory
   free(cp);
   free(c_pv);
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);
   free(ca_delay);
   free(wb_delay);
   free(c_t);
}

//------------------------------------------------------------------------------
// kconv_liver_jac
//------------------------------------------------------------------------------
// Calculates the time activity curves and sensitivity functions using the
// two-tissue compartment model with a dual-blood input function for the liver.
void kconv_liver_jac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
{
   int      i, j;
   double   vb, k1, k2, k3, k4, ka, fa, t_delay;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double   *c_pv, *cp, *c_a1, *c_a2, *c_f, *c_b, *c_t, *s_t;
   double   *ca_delay, *ca_grad_delay;
   int      num_time;
   int      num_par;

   // Determine the number of sensitive parameters
   num_par = 0;
   for (i = 0; i < 8; i++) {
      if (psens[i] == 1)
         ++num_par;
   }

   // Allocate memory for intermediate calculations
   num_time = (int) (scant[2*num_frm-1] / td);
   c_pv = (double*) malloc(sizeof(double) * num_time);
   cp   = (double*) malloc(sizeof(double) * num_time);
   c_a1 = (double*) malloc(sizeof(double) * num_time);
   c_a2 = (double*) malloc(sizeof(double) * num_time);
   c_f  = (double*) malloc(sizeof(double) * num_time);
   c_b  = (double*) malloc(sizeof(double) * num_time);

   ca_delay   = (double*) malloc(sizeof(double) * num_time);
   ca_grad_delay  = (double*) malloc(sizeof(double)*num_time);

   c_t  = (double*) malloc(sizeof(double) * num_time);
   s_t  = (double*) malloc(sizeof(double) * num_time * num_par);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {

      // Transform kinetic parameters into exponential parameters
      vb = p[0 + j * 8];
      k1 = p[1 + j * 8];
      k2 = p[2 + j * 8];
      k3 = p[3 + j * 8];
      k4 = p[4 + j * 8];
      ka = p[5 + j * 8];
      fa = p[6 + j * 8];
      t_delay = p[7 + j * 8];

      k234 = k2 + k3 + k4;
      d = sqrt(k234 * k234 - 4 * k2 * k4);
      a1 = (k234 - d) / 2;
      a2 = (k234 + d) / 2;

      // shift ca for time delay
      time_delay_tac(ca, num_time, t_delay, td, ca_delay);

      // Calculate dispersion
      tmp = ka + dk;
      kconv_exp(ka, tmp, ca_delay, num_time, td, c_pv);
      for (i = 0; i < num_time; i++) {
         cp[i] = (1 - fa) * c_pv[i] + fa * ca_delay[i];
      }

      // Calculate the exponential components
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2);

      // Calculate concentrations and sensitivity functions
      if (d == 0) d = 1e9;
      f1 = k1 / d * (k4 - a1);
      f2 = k1 / d * (a2 - k4);
      b1 = k1 / d * k3;
      b2 = -b1;
      for (i = 0; i < num_time; i++) {
         c_f[i] = f1 * c_a1[i] + f2 * c_a2[i];
         c_b[i] = b1 * c_a1[i] + b2 * c_a2[i];
         c_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * cp[i];
      }

      // Average the time activity curve over the frame duration
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);

      // Calculate sensitivity functions for each parameter
      if (psens[0] == 1) { // Sensitivity wrt vb
         for (i = 0; i < num_time; i++)
            s_t[i] = - (c_f[i] + c_b[i]) + cp[i];
         s_t += num_time;
      }
      if (psens[1] == 1 || psens[2] == 1) { 
         f1 = 1 / d * (k4 + k3 - a1);
         f2 = 1 / d * (a2 - k4 - k3);
      }
      if (psens[1] == 1) { // Sensitivity wrt k1
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[2] == 1 || psens[3] == 1) { 
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a2);
      }
      if (psens[2] == 1) { // Sensitivity wrt k2
         for (i = 0; i < num_time; i++)
            s_t[i] = -(1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[3] == 1 || psens[4] == 1) { 
         f1 = 1 / d * (a1 + a2 - k3 - k4);
         f2 = -f1;
      }
      if (psens[3] == 1) { // Sensitivity wrt k3
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[4] == 1) { // Sensitivity wrt k4
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a2);
         for (i = 0; i < num_time; i++)
            s_t[i] = -(1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[5] == 1 || psens[6] == 1) { 
         for (i = 0; i < num_time; i++) {
            cp[i] = ca_delay[i] - c_pv[i];
         }
         f1 = k1 / d * (k4 - a1);
         f2 = k1 / d * (a2 - k4);
      }
      if (psens[5] == 1) { // Sensitivity wrt ka
         tmp = ka + dk;
         kconv_exp(1.0, tmp, cp, num_time, td, c_t);
         for (i = 0; i < num_time; i++) {
            c_t[i] = (1.0 - fa) * c_t[i];
         }
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_t, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_t, num_time, td, c_a2);
         for (i = 0; i < num_time; i++) {
            s_t[i] = (1 - vb) * ((f1 + b1) * c_a1[i] + (f2 + b2) * c_a2[i]) + vb * c_t[i];
         }
         s_t += num_time;
      }
      if (psens[6] == 1) { // Sensitivity wrt fa
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, cp, num_time, td, c_a2);
         for (i = 0; i < num_time; i++) {
            s_t[i] = (1 - vb) * ((f1 + b1) * c_a1[i] + (f2 + b2) * c_a2[i]) + vb * cp[i];
         }
         s_t += num_time;
      }
      if (psens[7] == 1) { // wrt time delay
         time_delay_jac(ca_delay, num_time, td, ca_grad_delay);
         tmp = ka + dk;
         kconv_exp(ka, tmp, ca_grad_delay, num_time, td, c_pv);
         for (i = 0; i < num_time; i++) {
            cp[i] = (1 - fa) * c_pv[i] + fa * ca_grad_delay[i];
         }
         // Calculate the exponential components
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, cp, num_time, td, c_a2);

         // Calculate concentrations and sensitivity functions
         if (d == 0) d = 1e9;
         f1 = k1 / d * (k4 - a1);
         f2 = k1 / d * (a2 - k4);
         b1 = k1 / d * k3;
         b2 = -b1;

	      for (i=0; i<num_time; i++){ 
            c_f[i] = f1 * c_a1[i] + f2 * c_a2[i];
            c_b[i] = b1 * c_a1[i] + b2 * c_a2[i];
            s_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * cp[i];
            s_t[i] = - s_t[i];
         } 
	      s_t += num_time;
      }
      s_t -= num_time * num_par;

      // Average the sensitivity functions over the frame duration
      frame(scant, td, s_t, num_frm, num_par, st + j * num_frm * num_par);
   }

   // Free allocated memory
   free(cp);
   free(c_pv);
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);

   free(ca_delay);
   free(ca_grad_delay);

   free(c_t);
   free(s_t);
}

//------------------------------------------------------------------------------

