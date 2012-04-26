#include <math.h>
#include <stdio.h>

double erfsun_(double);

int main(){
  float a=0.12345;
  float b, c;

  b=erfsun_(a);
  c=erf(a);

  printf("%15.10f %15.10f %15.10f  \n",a,b,c);
  
  return 0;
}
