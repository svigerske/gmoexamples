#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "loadgms.h"

#include "gmomcc.h"
#include "gevmcc.h"
#include "GamsNLinstr.h"

static
RETURN dumpinstance(
   gmoHandle_t gmo,
   const char* filename
)
{
   FILE* file;
   char msg[GMS_SSSIZE];
   int rc;

   gmoOptFileSet(gmo, 1);
   gmoNameOptFileSet(gmo, "convertd.o99");
   file = fopen("convertd.o99", "w");
   fprintf(file, "gams %s\n", filename);
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

   return RETURN_OK;
}

static
RETURN copyGMO(
   gmoHandle_t  gmo,
   gmoHandle_t* gmocopy,
   int nextraconstants
   )
{
   gevHandle_t gev;
   char msg[GMS_SSSIZE];
   int i, j;
   int* colidx;
   double* jacval;
   int* nlflag;
   int nz, nlnz;
   int nconstants;
   double* constants;
   int ninstr;
   int* opcodes;
   int* fields;

   gev = gmoEnvironment(gmo);

   /* create empty GMO object */
   if( !gmoCreate(gmocopy, msg, sizeof(msg)) )
   {
      gevLogStatPChar(gev, "Could not create gmo object: ");
      gevLogStat(gev, msg);
      return RETURN_ERROR;
   }
   gmoRegisterEnvironment(*gmocopy, gev, msg);

   gmoObjStyleSet(gmo, gmoObjType_Var);

   /* attributes required to set before adding rows/cols */
   assert(gmoIndexBase(gmo) == 0);
   gmoIndexBaseSet(*gmocopy, 0);
   gmoObjStyleSet(*gmocopy, gmoObjType_Var);
   gmoPinfSet(*gmocopy, gmoPinf(gmo));
   gmoMinfSet(*gmocopy, gmoMinf(gmo));

   /* direction must be set before adding the matrix */
   gmoSenseSet(*gmocopy, gmoSense(gmo));

   /* model Type must be set before init data */
   gmoModelTypeSet(*gmocopy, gmoModelType(gmo));

   /* initialize model */
   gmoInitData(*gmocopy, gmoM(gmo) + 1, gmoN(gmo) + 1, 0);

   /* add all variables from gmo */
   for( j = 0; j < gmoN(gmo); ++j )
   {
      gmoAddCol(*gmocopy,
         gmoGetVarTypeOne(gmo, j),
         gmoGetVarLowerOne(gmo, j),
         gmoGetVarLOne(gmo, j),
         gmoGetVarUpperOne(gmo, j),
         gmoGetVarMOne(gmo, j),
         gmoGetVarStatOne(gmo, j),
         gmoGetVarSosSetOne(gmo, j),
         gmoGetVarPriorOne(gmo, j),
         gmoGetVarScaleOne(gmo, j),
         0, NULL, NULL, NULL
      );
   }

   gmoObjVarSet(*gmocopy, gmoObjVar(gmo));

   colidx = (int*)malloc(gmoN(gmo) * sizeof(int));
   jacval = (double*)malloc(gmoN(gmo) * sizeof(double));
   nlflag = (int*)malloc(gmoN(gmo) * sizeof(int));

   /* GMO stores one array with the constants that appear in nonlinear instructions: the constants pool
    * Once it has been set, it cannot be extended anymore (no API function).
    * Thus, we already set here the constants pool to the one from the original GMO followed by some extra space (nextraconstants)
    */
   nconstants = gmoNLConst(gmo) + nextraconstants;
   constants = (double*) malloc(nconstants * sizeof(double));
   memcpy(constants, gmoPPool(gmo), gmoNLConst(gmo) * sizeof(double));
   memset(constants + gmoNLConst(gmo), 0, nextraconstants * sizeof(double));

   ninstr = gmoNLCodeSizeMaxRow(gmo);
   opcodes = (int*) malloc(ninstr * sizeof(int));
   fields = (int*) malloc(ninstr * sizeof(int));

   for( i = 0; i < gmoM(gmo); ++i )
   {
      /* add equations from MINLP */

      gmoGetRowSparse(gmo, i, colidx, jacval, nlflag, &nz, &nlnz);
      gmoAddRow(*gmocopy,
         gmoGetEquTypeOne(gmo, i),
         0, /* gmoGetEquMatchOne(gmo, i), */
         gmoGetEquSlackOne(gmo, i),
         gmoGetEquScaleOne(gmo, i),
         gmoGetRhsOne(gmo, i),
         gmoGetEquMOne(gmo, i),
         gmoGetEquStatOne(gmo, i),
         nz, colidx, jacval, nlflag
      );

      /* copy nonlinear instructions
       * the first call to gmoDirtySetRowFNLInstr will set the constants pool in gmocopy
       */
      gmoDirtyGetRowFNLInstr(gmo, i, &ninstr, opcodes, fields);
      if( ninstr > 0 )
         gmoDirtySetRowFNLInstr(*gmocopy, i, ninstr, opcodes, fields, NULL, constants, nconstants);
   }

   /* more attributes before CompleteData (borrow from original model) */
   gmoPriorOptSet  (*gmocopy, gmoPriorOpt(gmo));
   gmoScaleOptSet  (*gmocopy, gmoScaleOpt(gmo));
   gmoDictionarySet(*gmocopy, gmoDictionary(gmo));
   gmoHaveBasisSet (*gmocopy, gmoHaveBasis(gmo));
   gmoNameXLib(gmo, msg); gmoNameXLibSet(*gmocopy, msg);
   gmoNameDict(gmo, msg); gmoNameDictSet(*gmocopy, msg);

   gmoCompleteData(*gmocopy, msg);


   free(opcodes);
   free(fields);
   free(constants);
   free(colidx);
   free(jacval);
   free(nlflag);

   return RETURN_OK;
}


/* add row  sum_i (x_i - \hat x_i) <= 10, where \hat x is the initial point */
static
RETURN addrow(
   gmoHandle_t gmo
)
{
   int* colidx;
   double* jacval;
   int* nlflag;
   int* opcodes;
   int* fields;
   int i;
   int ninstr;
   char msg[GMS_SSSIZE];


   /* Store current level values of all variables at the end of the constants pool. */
   gmoGetVarL(gmo, (double*)gmoPPool(gmo) + (gmoNLConst(gmo) - gmoN(gmo) + 1));

   /* setup instructions for sum_i sqr(x_i - \hat x_i)
    * this should be
    * pushv x_0
    * subi \hat x_0
    * callfunc sqr
    * pushv x_1
    * subi \hat x_1
    * callfunc sqr
    * add
    * pushv x_2
    * subi \hat x_2
    * callfunc sqr
    * add
    * ...
    *
    * setup also colidx, jacval, and nlflags for x_i's
    */

   ninstr = 4*gmoN(gmo)+1;
   opcodes = (int*) malloc(ninstr * sizeof(int));
   fields = (int*) malloc(ninstr * sizeof(int));

   colidx = (int*) malloc(gmoN(gmo) * sizeof(int));
   jacval = (double*) malloc(gmoN(gmo) * sizeof(double));
   nlflag = (int*) malloc(gmoN(gmo) * sizeof(int));

   opcodes[0] = nlHeader;
   fields[0] = 0;
   ninstr = 1;

   for( i = 0; i < gmoN(gmo); ++i )
   {
      opcodes[ninstr] = nlPushV;
      fields[ninstr] = gmoGetjModel(gmo, i) + 1;  /* use index in "model-space" (GMO shows "solver-space" to user) and +1 for 1-based indexing (??) */
      ++ninstr;

      opcodes[ninstr] = nlSubI;
      fields[ninstr] = gmoNLConst(gmo) - gmoN(gmo) + i + 2;  /* +1 as GMO looks into constants pool in 1-based indexing */
      ++ninstr;

      opcodes[ninstr] = nlCallArg1;
      fields[ninstr] = fnsqr;
      ++ninstr;

      if( i > 0 )
      {
         opcodes[ninstr] = nlAdd;
         fields[ninstr] = 0;
         ++ninstr;
      }

      colidx[i] = i;   /* these are in "solver-space" again */
      jacval[i] = 1.0; /* arbitrary, as nonlinear */
      nlflag[i] = 1;   /* indicate that nonlinear */
   }

   opcodes[ninstr] = nlStore;
   fields[ninstr] = 0;
   ++ninstr;
   assert(ninstr == 4*gmoN(gmo)+1);

   /* add row */
   gmoAddRow(gmo, gmoequ_L, 0, 0.0, 0.0, 10.0 /*rhs*/, 0.0, gmoBstat_Upper, gmoN(gmo), colidx, jacval, nlflag);
   gmoDirtySetRowFNLInstr(gmo, gmoM(gmo)-1, ninstr, opcodes, fields, NULL, NULL, 0);

   gmoCompleteData(gmo, msg);

   free(opcodes);
   free(fields);
   free(colidx);
   free(jacval);
   free(nlflag);

   return RETURN_OK;
}


int main(int argc, char** argv)
{
   gmoHandle_t gmo;
   gmoHandle_t gmocopy;
   gevHandle_t gev;

   if( argc < 2 )
   {
      printf("usage: %s <.gms file>\n", argv[0]);
      return EXIT_FAILURE;
   }

   CHECK( loadGMS(&gmo, &gev, argv[1]) );

   CHECK( copyGMO(gmo, &gmocopy, gmoN(gmo)) );

   CHECK( dumpinstance(gmocopy, "before.gms") );

   CHECK( addrow(gmocopy) );

   CHECK( dumpinstance(gmocopy, "after.gms") );

   system("diff before.gms after.gms");

   gmoFree(&gmocopy);

   freeGMS(&gmo, &gev);

   return EXIT_SUCCESS;
}
