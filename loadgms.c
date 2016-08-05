#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "gmomcc.h"
#include "gevmcc.h"
#include "assert.h"

#include "loadgms.h"

#define GAMSLOGOPTION 0

RETURN loadGMS(
   struct gmoRec** gmo,
   struct gevRec** gev,
   const char* gmsfile
)
{
   char gamscall[1024];
   char buffer[GMS_SSSIZE];
   int rc;
   
   FILE* convertdopt;
   
   /* create temporary directory */
   mkdir("loadgms.tmp", S_IRWXU);
   
   /* create empty convertd options file */
   convertdopt = fopen("loadgms.tmp/convertd.opt", "w");
   if( convertdopt == NULL )
   {
      fprintf(stderr, "Could not create convertd options file.\n");
      return RETURN_ERROR;
   }
   fputs(" ", convertdopt);
   fclose(convertdopt);
   
   /* call GAMS with convertd solver to get compiled model instance in temporary directory */
   snprintf(gamscall, sizeof(gamscall), GAMSDIR "/gams %s LP=CONVERTD RMIP=CONVERTD QCP=CONVERTD RMIQCP=CONVERTD NLP=CONVERTD DNLP=CONVERTD RMINLP=CONVERTD CNS=CONVERTD MIP=CONVERTD MIQCP=CONVERTD MINLP=CONVERTD MCP=CONVERTD MPEC=CONVERTD RMPEC=CONVERTD SCRDIR=loadgms.tmp output=loadgms.tmp/listing optdir=loadgms.tmp optfile=1 pf4=0 solprint=0 limcol=0 limrow=0 pc=2 lo=%d", gmsfile, GAMSLOGOPTION);
   /* printf(gamscall); fflush(stdout); */
   rc = system(gamscall);
   if( rc != 0 )
   {
      fprintf(stderr, "GAMS call returned with code %d\n", rc);
      return RETURN_ERROR;
   }

   /* initialize GMO and GEV libraries */
   if( !gmoCreateDD(gmo, GAMSDIR, buffer, sizeof(buffer)) || !gevCreateDD(gev, GAMSDIR, buffer, sizeof(buffer)) )
   {
      fprintf(stderr, buffer);
      fprintf(stderr, "\n");
      return 1;
   }

   /* load control file */
   if( gevInitEnvironmentLegacy(*gev, "loadgms.tmp/gamscntr.dat") )
   {
      fprintf(stderr, "Could not load control file loadgms.tmp/gamscntr.dat\n");
      gmoFree(gmo);
      gevFree(gev);
      return 1;
   }

   if( gmoRegisterEnvironment(*gmo, *gev, buffer) )
   {
      fprintf(stderr, "Error registering GAMS Environment: %s\n", buffer);
      gmoFree(gmo);
      gevFree(gev);
      return 1;
   }

   if( gmoLoadDataLegacy(*gmo, buffer) )
   {
      fprintf(stderr, "Could not load model data.\n");
      gmoFree(gmo);
      gevFree(gev);
      return 1;
   }

   /* reformulate objective variable out of model, if possible */
   gmoObjStyleSet(*gmo, gmoObjType_Fun);
   gmoObjReformSet(*gmo, 1);

   gmoIndexBaseSet(*gmo, 0);

   return RETURN_OK;
}

void freeGMS(
   struct gmoRec** gmo,
   struct gevRec** gev
)
{
   gmoFree(gmo);
   gmo = NULL;
   
   gevFree(gev);
   gev = NULL;
   
   gmoLibraryUnload();
   gevLibraryUnload();
   
   /* remove temporary directory content (should have only files) and directory itself) */
   system("rm loadgms.tmp/* && rmdir loadgms.tmp");
}
