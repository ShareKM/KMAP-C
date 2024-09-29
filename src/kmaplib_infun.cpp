/*
 * kmaplib_infun.cpp
 * 
 * This file contains the functions that are used for processing the input function.
 */

#include "kmaplib.h"
#include <cmath>
#include <cstdlib>

//------------------------------------------------------------------------------
// time_delay_tac
//------------------------------------------------------------------------------
// compute time-delayed TAC curve
void time_delay_tac(double* tac, int tac_size, double delay_time, double td, double *out) {
/* 
   tac: time activity curve
   tac_size: number of frames
   delay_time: time delay (in seconds)
   td: time step size
   out: output tac

   Created by Yiran Wang, modified by Yansong Zhu @ UC Davis
*/
   double shiftAmount = delay_time/td;
   double fraction = shiftAmount - floor(shiftAmount);
   int shift = floor(shiftAmount);
   // Handling right shift
   if (shift >= 0) {
      int x1;
      for(int i = 0; i < min(shift, tac_size); i++) { //padding right side with 0
         x1 = i;
         out[x1] = 0;
      }
      if (shift < tac_size) { //compute values of the shifted left boundary point with linear interpolation
         out[shift] = (1-fraction) * tac[0];   
      }

      int z1,z2,z3; 
      for(int i = shift; i < tac_size-1; i++){ //compute values of the tac curve for other sample points with linear interpolation
         z1=i+1; z2=i-shift; z3=i-shift+1;
         out[z1] = (fraction) * tac[z2] + (1-fraction) * tac[z3]; 
      }
   } 
   // Handling left shift
   else {
      shift = -shift; // make shift positive for easier handling
      shift = shift - 1;

      int a1,a2,a3;
      for(int i = 0; i < tac_size - shift - 1; i++) { //compute values of the tac curve for other sample points with linear interpolation
         a1=i; a2=i+shift; a3=i+shift+1;
         out[a1] = (fraction) * tac[a2] + (1-fraction) * tac[a3];
      }
      if (shift < tac_size) { //compute values of the shifted right boundary point with linear interpolation
         int b1, b2, b3; b1=tac_size-1-shift; b2=tac_size-1; b3=tac_size-1;
         out[b1] = (fraction) * tac[b2] + (1-fraction) * tac[b3]; 
      }
        
      int c1,c2;
      for(int i = tac_size - shift; i < tac_size; i++) {
         c1=i; c2=tac_size-1;
         out[c1] = tac[c2]; // padding new points with the tail value of the original data
      }
   }
}

//------------------------------------------------------------------------------
// time_delay_jac
//------------------------------------------------------------------------------
// compute gradient for time delay correction
void time_delay_jac(double *tac, int tac_size, double delay_time, double td, double *out){
/*
   tac: time activity curve
   tac_size: number of frames
   delay_time: time delay (in seconds)
   td: time step size
   out: output gradient vector

   Created by Yiran Wang in 2023, modified by Yansong Zhu in 2024 @ UC Davis
*/
int n = floor(delay_time/td);
for(int m=0; m<tac_size; m++){
   int i = m - (n + 1);
   if ((i >= 0) && (i <= (tac_size - 2))){
      out[m] = (tac[i + 1] - tac[i]) / td;
   } 
   else{
      out[m] = 0.0;    
   }
}
}
