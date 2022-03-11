#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#include "gmomcc.h"
#include "gevmcc.h"
#include "assert.h"

#include "loadgms.h"

/* set to 3 to see gams log */
#define GAMSLOGOPTION 0

RETURN loadGMS(
   struct gmoRec** gmo,
   struct gevRec** gev,
   const char*     gmsfile
)
{
   char scrdir[] = "/tmp/loadgmsXXXXXX";
   char filename[35];
   char gamscall[1024];
   char buffer[GMS_SSSIZE];
   FILE* convertopt;
   int rc;
   
   /* create temporary directory */
   mkdtemp(scrdir);
   
   /* create convert options file */
   sprintf(filename, "%s/convert.opt", scrdir);
   convertopt = fopen(filename, "w");
   if( convertopt == NULL )
   {
      fprintf(stderr, "Could not create convert options file %s.\n", filename);
      return RETURN_ERROR;
   }
   /* with new convert, something needs to be created to avoid writing gams.gms and dict.txt into working dir */
   fprintf(convertopt, "dict %s/dict\n", scrdir);
   fclose(convertopt);
   
   /* call GAMS with convert solver to get compiled model instance in temporary directory */
   snprintf(gamscall, sizeof(gamscall), "\"" GAMSDIR "/gams\" %s SOLVER=CONVERT SCRDIR=%s output=%s/listing optdir=%s optfile=1 pf4=0 solprint=0 limcol=0 limrow=0 pc=2 lo=%d", gmsfile, scrdir, scrdir, scrdir, GAMSLOGOPTION);
   rc = system(gamscall);
   if( rc != 0 )
   {
      fprintf(stderr, "GAMS call returned with code %d\n", rc);
      return RETURN_ERROR;
   }

   /* initialize GMO and GEV libraries */
   if( !gmoCreateDD(gmo, GAMSDIR, buffer, sizeof(buffer)) || !gevCreateDD(gev, GAMSDIR, buffer, sizeof(buffer)) )
   {
      fputs(buffer, stderr);
      fputs("\n", stderr);
      return 1;
   }

   /* load control file */
   sprintf(filename, "%s/gamscntr.dat", scrdir);
   if( gevInitEnvironmentLegacy(*gev, filename) )
   {
      fprintf(stderr, "Could not load control file %s\n", filename);
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
   char scrdirname[20];
   char filename[275];
   DIR* scrdir;
   struct dirent* scrdirentry;

   gevGetStrOpt(*gev, gevNameScrDir, scrdirname);

   gmoFree(gmo);
   *gmo = NULL;
   
   gevFree(gev);
   *gev = NULL;
   
   gmoLibraryUnload();
   gevLibraryUnload();
   
   /* remove scratch dir content: should be only files */
   scrdir = opendir(scrdirname);
   if( scrdir != NULL )
   {
      while( (scrdirentry = readdir(scrdir)) != NULL )
      {
         if( scrdirentry->d_type != DT_DIR )  /* skip . and .. */
         {
            sprintf(filename, "%s%s", scrdirname, scrdirentry->d_name);
            /* printf("Removing %s\n", filename); */
            if( remove(filename) != 0 )
               fprintf(stderr, "Could not remove temporary file %s\n", filename);
         }
      }
   }
   closedir(scrdir);

   /* remove scratch dir: should be empty now */
   if( remove(scrdirname) != 0 )
      fprintf(stderr, "Could not remove temporary directory %s\n", scrdirname);
}
