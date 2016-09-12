#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "loadgms.h"

#include "gmomcc.h"
#include "gevmcc.h"
#include "GamsNLinstr.h"

/** Writes model instance into a .gms file of given name.
 *
 * The function calls GAMS/ConvertD to write the instance.
 */
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

/** Creates a copy of a GMO object with the possibility to extended the constants pool.
 *
 * Matchings (MCP) are not copied so far.
 */
static
RETURN copyGMO(
   gmoHandle_t  gmo,             /**< GMO to copy */
   gmoHandle_t* gmocopy,         /**< buffer to store pointer to new GMO */
   int          nextraconstants  /**< the number of extra constants for which space in the constants pool of gmocopy should be reserved */
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

   /* only know how to copy if objective variable has not been reformulated out */
   gmoObjStyleSet(gmo, gmoObjType_Var);

   /* attributes required to set before adding rows/cols */
   assert(gmoIndexBase(gmo) == 0);
   gmoIndexBaseSet(*gmocopy, 0);
   gmoObjStyleSet(*gmocopy, gmoObjType_Var);
   gmoPinfSet(*gmocopy, gmoPinf(gmo));
   gmoMinfSet(*gmocopy, gmoMinf(gmo));

   /* direction must be set before adding the matrix */
   gmoSenseSet(*gmocopy, gmoSense(gmo));

   /* model type must be set before init data */
   gmoModelTypeSet(*gmocopy, gmoModelType(gmo));

   /* initialize model */
   gmoInitData(*gmocopy, gmoM(gmo), gmoN(gmo), 0);

   /* add all variables from gmo */
   for( j = 0; j < gmoN(gmo); ++j )
   {
      gmoAddCol(*gmocopy,
         gmoGetVarTypeOne(gmo, j),
         gmoGetVarLowerOne(gmo, j),
         gmoGetVarLOne(gmo, j),
         gmoGetVarUpperOne(gmo, j),
         gmoGetVarMOne(gmo, j),
         gmoGetVarStatOne(gmo, j) != gmoBstat_Basic,
         gmoGetVarSosSetOne(gmo, j),
         gmoGetVarPriorOne(gmo, j),
         gmoGetVarScaleOne(gmo, j),
         0, NULL, NULL, NULL
      );
   }

   /* say which one is the objective variable */
   gmoObjVarSet(*gmocopy, gmoObjVar(gmo));

   /* space for coefficients and nonlinear instructions in each row */
   colidx = (int*)malloc(gmoN(gmo) * sizeof(int));
   jacval = (double*)malloc(gmoN(gmo) * sizeof(double));
   nlflag = (int*)malloc(gmoN(gmo) * sizeof(int));
   opcodes = (int*) malloc(gmoNLCodeSizeMaxRow(gmo) * sizeof(int));
   fields = (int*) malloc(gmoNLCodeSizeMaxRow(gmo) * sizeof(int));

   /* GMO stores one array with the constants that appear in nonlinear instructions: the constants pool
    * Once it has been set, it cannot be extended anymore (no API function).
    * Thus, we already set here the constants pool to the one from the original GMO followed by some extra space (nextraconstants)
    */
   nconstants = gmoNLConst(gmo) + nextraconstants;
   constants = (double*) malloc(nconstants * sizeof(double));
   memcpy(constants, gmoPPool(gmo), gmoNLConst(gmo) * sizeof(double));
   memset(constants + gmoNLConst(gmo), 0, nextraconstants * sizeof(double));

   /* add all rows from gmo */
   for( i = 0; i < gmoM(gmo); ++i )
   {
      /* copy coefficients */
      gmoGetRowSparse(gmo, i, colidx, jacval, nlflag, &nz, &nlnz);
      gmoAddRow(*gmocopy,
         gmoGetEquTypeOne(gmo, i),
         0, /* gmoGetEquMatchOne(gmo, i), */
         gmoGetEquSlackOne(gmo, i),
         gmoGetEquScaleOne(gmo, i),
         gmoGetRhsOne(gmo, i),
         gmoGetEquMOne(gmo, i),
         gmoGetEquStatOne(gmo, i) != gmoBstat_Basic,
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


/* add row  sum_{i!=objvar} sqr(x_i - \hat x_i) <= 10, where \hat x is the initial point */
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

   /* we will need 4*(nvars-1) + 1 instructions */
   ninstr = 4*(gmoN(gmo)-1)+1;
   opcodes = (int*) malloc(ninstr * sizeof(int));
   fields = (int*) malloc(ninstr * sizeof(int));

   colidx = (int*) malloc((gmoN(gmo)-1) * sizeof(int));
   jacval = (double*) malloc((gmoN(gmo)-1) * sizeof(double));
   nlflag = (int*) malloc((gmoN(gmo)-1) * sizeof(int));

   opcodes[0] = nlHeader;
   /* fields[0] will be set to the number of instructions below */;
   ninstr = 1;

   for( i = 0; i < gmoN(gmo); ++i )
   {
      int ii;
      // skip objective variable
      // including it here will prevent restoring the objective function form the objective equation,
      // which GMO assumes that it can do even after adding additional rows including the objective variable
      // however, calls like gmoGetMatrixRow get confused on the NZ count
      if( i == gmoObjVar(gmo) )
         continue;

      opcodes[ninstr] = nlPushV;
      fields[ninstr] = gmoGetjModel(gmo, i) + 1;  /* use index in "model-space" (GMO shows "solver-space" to user) and +1 for 1-based indexing (??) */
      ++ninstr;

      opcodes[ninstr] = nlSubI;
      fields[ninstr] = gmoNLConst(gmo) - gmoN(gmo) + i + 2;  /* +1 as GMO looks into constants pool in 1-based indexing, and another +1 */
      ++ninstr;

      opcodes[ninstr] = nlCallArg1;
      fields[ninstr] = fnsqr;
      ++ninstr;

      if( ninstr > 4 )
      {
         opcodes[ninstr] = nlAdd;
         fields[ninstr] = 0;
         ++ninstr;
      }

      ii = (i > gmoObjVar(gmo) ? i-1 : i);
      colidx[ii] = i;   /* these are in "solver-space" again */
      jacval[ii] = 1.0; /* arbitrary, as nonlinear */
      nlflag[ii] = 1;   /* indicate that nonlinear */
   }

   opcodes[ninstr] = nlStore;
   fields[ninstr] = gmoM(gmo);
   ++ninstr;
   assert(ninstr == 4*(gmoN(gmo)-1)+1);

   fields[0] = ninstr;

   /* add row */
   gmoAddRow(gmo, gmoequ_L, 0, 0.0, 0.0, 10.0 /*rhs*/, 0.0, 0, gmoN(gmo)-1, colidx, jacval, nlflag);
   gmoDirtySetRowFNLInstr(gmo, gmoM(gmo)-1, ninstr, opcodes, fields, NULL, NULL, 0);

   free(opcodes);
   free(fields);
   free(colidx);
   free(jacval);
   free(nlflag);

   /* call completeData to update some additional counts
    * without this, equnextNL does not get updated in GMO, which will
    * prevent including the added row into the Lagrangian Hessian
    */
   if( gmoCompleteData(gmo, msg) )
   {
      gevLogStatPChar(gmoEnvironment(gmo), "Failure from gmoCompleteData: ");
      gevLogStat(gmoEnvironment(gmo), msg);
      return RETURN_ERROR;
   }

   return RETURN_OK;
}

static
void printHessian(
   gmoHandle_t gmo
)
{
   int* rowindex;
   int* colindex;
   int dim;
   int nz;
   int i;

   dim = gmoHessLagDim(gmo);
   nz = gmoHessLagNz(gmo);
   rowindex = (int*)malloc(nz * sizeof(int));
   colindex = (int*)malloc(nz * sizeof(int));

   gmoHessLagStruct(gmo, rowindex, colindex);

   printf("Hessian (%d nonzeros, dim %d):", nz, dim);
   for( i = 0; i < nz; ++i )
      printf(" (%d,%d)", rowindex[i], colindex[i]);
   printf("\n");
   printf("NLM: %d\n", gmoNLM(gmo));

   free(rowindex);
   free(colindex);
}


int main(int argc, char** argv)
{
   gmoHandle_t gmo;
   gmoHandle_t gmocopy;
   gevHandle_t gev;
   char msg[255];

   if( argc < 2 )
   {
      printf("usage: %s <.gms file>\n", argv[0]);
      return EXIT_FAILURE;
   }

   /* load .gms file into GMO */
   CHECK( loadGMS(&gmo, &gev, argv[1]) );

   /* Now we would like to add a nonlinear row to gmo.
    * However, there are (at least) two stumbling blocks that GMO puts in our way:
    * 1. The constants pool cannot be extended in GMO, as far as I see.
    *    gmoSetNLObject() gets close to it, but we don't know how to build the nlobject and nlpool objects.
    * 2. GMO holds one array to store all nonlinear instructions.
    *    gmoDirtySetRowFNLInstr does not add new nonlinear instructions to the end of this array, but
    *    where it stopped adding instructions on the previous call (variable "dirtydirty", initialized to 0).
    *    Thus, we can only add new nonlinear instructions, if we have also added all previous instructions.
    *    Further, we cannot reset the instructions of existing nonlinear rows to get the "dirtydirty" correct.
    *
    * Therefore, we create a copy of the instance in a new GMO object.
    * When doing this, we already reserve extra space in the constants pool.
    * Further, setting up the GMO object by row-by-row will update the "dirtydirty" counter so that it will
    * finally point to the end of the nonlinear instructions array.
    */
   CHECK( copyGMO(gmo, &gmocopy, gmoN(gmo)) );

   CHECK( dumpinstance(gmocopy, "before.gms") );

/*
   gevCallSolver(gmoEnvironment(gmocopy), gmocopy, "", "ipopt", gevSolveLinkLoadLibrary,
      gevSolverSameStreams, "", "",
      1000.0, ITERLIM_INFINITY, 0,
      0.0, 0.0, NULL, msg);
   printHessian(gmocopy);
*/

   CHECK( addrow(gmocopy) );

   CHECK( dumpinstance(gmocopy, "after.gms") );

   system("diff before.gms after.gms");

/*
   gevCallSolver(gmoEnvironment(gmocopy), gmocopy, "", "ipopt", gevSolveLinkLoadLibrary,
      gevSolverSameStreams, "", "",
      1000.0, ITERLIM_INFINITY, 0,
      0.0, 0.0, NULL, msg);
   printHessian(gmocopy);
*/

   gmoFree(&gmocopy);

   freeGMS(&gmo, &gev);

   return EXIT_SUCCESS;
}
