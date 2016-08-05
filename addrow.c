#include <stdio.h>
#include <stdlib.h>

#include "loadgms.h"

#include "gmomcc.h"
#include "gevmcc.h"

int main(int argc, char** argv)
{
   gmoHandle_t gmo;
   gevHandle_t gev;
   char buffer[1024];

   if( argc < 2 )
   {
      printf("usage: %s <.gms file>\n", argv[0]);
      return EXIT_FAILURE;
   }

   CHECK( loadGMS(&gmo, &gev, argv[1]) );

   printf("N vars: %d  rows: %d\n", gmoN(gmo), gmoM(gmo));



   return EXIT_SUCCESS;
}
