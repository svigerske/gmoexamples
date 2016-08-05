#ifndef DEF_H
#define DEF_H

typedef enum
{
    RETURN_OK = EXIT_SUCCESS,
    RETURN_ERROR = EXIT_FAILURE
} RETURN;


#define CHECK( x ) \
   do \
   { \
      RETURN _retcode = (x); \
      if( _retcode != RETURN_OK ) \
      { \
         fprintf(stderr, __FILE__ ":%d Error %d from call " #x "\n", __LINE__, _retcode); \
         return _retcode; \
      } \
   } while( 0 )

struct gmoRec;
struct gevRec;

#endif
