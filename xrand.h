/*

extern void init_random(int seed,int seed2);

   inizializza con init_random(0,mype) per avere una sequenza diversa
   ogni secondo e diversa per i diversi processori

extern double Xrandom(void);

   restituisce un numero random tra 0 e 1

*/

#include <time.h>

static int A[4];
static int B[4]={6698,7535,26792,30140};
                          /* 246913578*(2^32+1)=1060485742695258666 */

#define M0   17917        /* 13^13=302875106592253 */
#define M1   13895
#define M2   19930
#define M3   8

#define PW2  32767         /* 2^(15)-1 */

#define rpw2a  9.31322574615478516e-10     /* 2^(-30) */
#define rpw2b  1.11022302462515654e-16     /* 2^(-53) */
#define rpw2c  2.77555756156289135e-17     /* 2^(-55) */

inline double Xrandom(void)
  {
   A[0]=B[0]*M0;
   A[1]=(A[0]>>15)+B[1]*M0+B[0]*M1;
   A[2]=(A[1]>>15)+B[2]*M0+B[1]*M1+B[0]*M2;
   A[3]=(A[2]>>15)+B[3]*M0+B[2]*M1+B[1]*M2+B[0]*M3;
   B[0]=A[0]&PW2;
   B[1]=A[1]&PW2;
   B[2]=A[2]&PW2;
   B[3]=A[3]&PW2;
   return rpw2a*((B[3]<<15)+B[2])+rpw2b*((B[1]<<8)+(B[0]>>7))+rpw2c;
  }

void init_random(int seed,int seed2)
  {
   if (seed==0) seed=(int)time(NULL);
   B[0]=(2+(seed<<2))&PW2;
   B[1]=(seed>>13^seed2)&PW2;
   B[2]=(seed>>28^seed2>>15)&PW2;
   B[3]=seed2>>30&PW2;
   Xrandom();
  }

// gaussian number generator
bool hai_stored_gaussian = false;
double stored_gaussian;
double old_var = -1;

inline double Xrandom_gaussian(){
  // var2 = sigma^2
  float var2 = 1;
  if(hai_stored_gaussian && (old_var == var2)){
    hai_stored_gaussian = false;
    return stored_gaussian;
  }else{
//    old_var = var2;
    double v1,v2,r;
    do{
      v1 = 2*Xrandom() - 1;
      v2 = 2*Xrandom() - 1;
      r = v1*v1+v2*v2;
    }while( (r >= 1) || (r == 0));
    double fac = sqrt(-2*log(r)/r);
//    double fac = sqrt(-2*log(r)*var2/r);
    stored_gaussian = v1*fac;
    hai_stored_gaussian = true;
    return (v2*fac);
  }
}

inline double Xrandom_exponential(double lambda){
// distribuzione exp(-lambda*x)	
  double x;
  do{
    x = Xrandom();	
  }while(x == 0);
  return -log(x)*lambda;
}


