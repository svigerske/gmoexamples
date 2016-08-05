#ifndef LOADGMS_H
#define LOADGMS_H

#include "def.h"

extern
RETURN loadGMS(
   struct gmoRec** gmo,
   struct gevRec** gev,
   const char* gmsfile
);


extern
void freeGMS(
   struct gmoRec** gmo,
   struct gevRec** gev
);

#endif
