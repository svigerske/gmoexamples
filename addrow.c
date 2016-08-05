#include <stdio.h>
#include <stdlib.h>

#include "loadgms.h"

#include "gmomcc.h"
#include "gevmcc.h"

RETURN dumpinstance(
   gmoHandle_t gmo
)
{
   FILE* file;
   char msg[GMS_SSSIZE];
   int rc;

   gmoOptFileSet(gmo, 1);
   gmoNameOptFileSet(gmo, "convertd.o99");
   file = fopen("loadgms.tmp/convertd.o99", "w");
   fprintf(file, "gams dump.gms\n");
   fclose(file);

   rc = gevCallSolver(gmoEnvironment(gmo), gmo, "", "convertd", gevSolveLinkLoadLibrary,
      gevSolverSameStreams, "", "",
      1000.0, ITERLIM_INFINITY, 0,
      0.0, 0.0, NULL, msg);

   if( rc != 0 )
   {
      puts(msg);
      return RETURN_ERROR;
   }

   system("cat dump.gms");

   return RETURN_OK;
}

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

   CHECK( dumpinstance(gmo) );





   return EXIT_SUCCESS;
}
