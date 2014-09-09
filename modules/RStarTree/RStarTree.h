/* -----  RStarTree.h  ----- */
#ifndef __RStarTree_h
#define __RStarTree_h



/**                             R*-tree
                                =======                                 **/


/**    Implementation:  Norbert Beckmann
              Version:  R.2.0
                 Date:  6/93                                            **/


/**                     Praktische Informatik,
             Universitaet Bremen, D-2800 Bremen 33, Germany             **/



/* ---------------------- operating system version --------------------- */
/*
#ifndef SVR4
#  define SVR4
#endif
*/
#include <stddef.h>
#include <sys/file.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <sys/stat.h>

#include <fcntl.h>


/* padding on 32 bits systems (to match automatic 64 bits padding) */

#if defined (__LP64__) || defined (__64BIT__) || defined (_LP64) || (__WORDSIZE == 64)
  #define PADDING_32_BITS
#else
  #define PADDING_32_BITS int padding
#endif

/* ----------------------------- constants ----------------------------- */

#define byte unsigned char
#define boolean int
#define FALSE 0
#define TRUE 1

#define int int /**V**/ /* large data sets require a 32 bit int */

#define MaxNameLength 160 /* including less than 10 bytes for a suffix */

#define MaxNumbOfEntriesPerPage 512 /**V**/

#define NumbOfDim 2 /**V**/


/* ------------------------------- types ------------------------------- */

typedef int File;


typedef float  typatomkey; /**V**/

        /* typatomkey may be of any type as far as the standard comparisons
           apply */


typedef struct {
               typatomkey  l, h;
               } typinterval;


typedef typinterval  typrect[NumbOfDim];

        /* A typrect is the key type an R*-tree handles. The smallest
           entity which may be used as a key is an interval,
           i.e. a typrect[1]. */


typedef struct {
               float height;
               } typinfo, *refinfo;

        /* A typinfo is a struct which may contain arbitrary information
           associated with a data record.
           RESTRICTION: sizeof(typinfo) >= sizeof(int) must hold! */

typedef struct {
  double m01, m02, m03;
  double m11, m13;
  double m22, m23, m33;
  double m44, m55, m66, m77;
  double m67, m76;
  double H0, H1, H2, H3, H4;
  double H5, H6;
  float Hmin, Hmax;
  int n;
  PADDING_32_BITS;
} typdirinfo;

/* A typdirinfo is a struct which may contain arbitrary information
   associated with a directory record. */

typedef int (* Check) (typrect rect, void * data, int depth);

/* ------------------------- private includes -------------------------- */

#include "RSTBase.h"


/* ----------------------------- R*-tree ------------------------------- */

typedef rstree  *RSTREE; /* R*-tree identifier */


/* ------------------------- procedure types --------------------------- */

typedef boolean (*DirQueryProc) (RSTREE   /* rst */,
                                 typrect  /* dirrect */,
                                 typrect  /* queryrect1 */,
                                 typrect  /* queryrect2 */);

        /* rst:
           contains the R*-tree identifier as passed to the procedures
           RegionQuery, ExistsRegion and RegionCount respectively.
           
           dirrect:
           contains a multidimensional directory rectangle.
           
           queryrect1, queryrect2:
           contain queryrect1 and queryrect2 as passed to the procedures
           RegionQuery, ExistsRegion and RegionCount respectively.
           
           This type is also used in the JoinXX procedures to perform the
           join. Then the parameters will have a slightly different
           contents, see the JoinXX procedures. */


typedef boolean (*DataQueryProc) (RSTREE   /* rst */,
                                  typrect  /* datarect */,
                                  typrect  /* queryrect1 */,
                                  typrect  /* queryrect2 */);

        /* rst:
           contains the R*-tree identifier as passed to the procedures
           RegionQuery, ExistsRegion and RegionCount respectively.
           
           datarect:
           contains a multidimensional data rectangle.
           
           queryrect1, queryrect2:
           contain queryrect1 and queryrect2 as passed to the procedures
           RegionQuery, ExistsRegion and RegionCount respectively.
           
           This type is also used in the JoinXX procedures to perform the
           join. Then the parameters will have a slightly different
           contents, see the JoinXX procedures. */


typedef void (*QueryManageProc) (RSTREE    /* rst */,
                                 typrect   /* rectangle */,
                                 refinfo   /* info */,
                                 void*     /* pointer */,
                                 boolean*  /* modify */,
                                 boolean*  /* finish */);

        /* rst:
           contains the R*-tree identifier as passed to the procedure
           RegionQuery.
           
           rectangle:
           contains the multidimensional data rectangle currently found.
           
           info:
           points to the concerning information part.
           
           pointer:
           is an arbitrary pointer as passed to the procedure RegionQuery.
           It may be used as a pointer to a buffer structure.
           
           modify:
           points to a flag labeling the node to be written back.
           BEWARE:
           record modified  |  modify  |  result
           --------------------------------------------------------
                 no         |  FALSE   |  record is not modified
                 no         |  TRUE    |  undefined
                 yes        |  FALSE   |  undefined
                 yes        |  TRUE    |  record is modified
           (See also RegionQuery)
           
           finish:
           points to a flag labeling the query to be finished when the
           procedure is left. */


typedef void (*JoinManageProc) (RSTREE    /* rst1 */,
                                RSTREE    /* rst2 */,
                                typrect   /* rectangle1 */,
                                typrect   /* rectangle2 */,
                                refinfo   /* info1 */ ,
                                refinfo   /* info2 */,
                                void*     /* pointer1 */,
                                void*     /* pointer2 */,
                                boolean*  /* finish */);

        /* rst1, rst2:
           contain the R*-tree identifiers as passed to the JoinXX
           procedures.
           
           rectangle1, rectangle2:
           contain the two rectangles currently found by the join.
           
           info1, info2:
           point to the concerning information parts.
           
           pointer1, pointer2:
           are arbitrary pointers as passed to the JoinXX procedures.
           They may be used as pointers to buffer structures.
           
           finish:
           points to a flag labeling the join to be finished when the
           procedure is left. */


/* ---------------------- procedure declarations ----------------------- */

/*
              Almost all procedures return a boolean result.
     If a procedure returns FALSE an error has occured. The implementation
     attempts to detect errors before any update operations are performed.
*/


void     NoRSTree(RSTREE *rst);

         /* Initializes an R*-tree identifier (at least sets *rst to NULL).
            Each procedure which requires an R*-tree identifier checks the
            value of this identifier at its entry.
            OpenRST for example demands a NULL identifier while CloseRST
            demands a non NULL identifier. */


boolean  CreateRST(const char     *name,
                   int      pagesize,
                   boolean  unique);

         /* CreateRST creates an R*-tree on secondary memory.
            To work on it, it has to be opened by the procedure OpenRST.
            
            name:
            is the main filename under which the created R*-tree will be
            stored, additionally a few files named filename.suffix with
            different suffixes will be created. name is not fixed in the
            internal parameter table, thus after renaming the files, the
            R*-tree may be opened under another name.
            
            pagesize:
            is the size in bytes, a page (directory or data) will occupy.
            
            unique:
            if unique is set TRUE the procedure InsertRecord will not store
            more than one record with the same rectangle (key) in the
            R*-tree (rectangles will be real keys).
            The unique flag may be reset with the procedure SetUnique
            without further internal tests. */


boolean  RemoveRST(const char *name);

         /* RemoveRST removes all files corresponding to an R*-tree. */


boolean  OpenRST(RSTREE        *rst,
                 const char    *name,
		 const char    *mode);

         /* OpenRST opens the R*-tree named name. */


boolean  CloseRST(RSTREE *rst);

         /* CloseRST closes the R*-tree referenced by rst. */


boolean  SetUnique(RSTREE   rst,
                   boolean  mode);

         /* The unique state, defined in procedure CreateRST may be reset
            by this procedure (see also CreateRST). The unique flag is set
            without internal checks (even to TRUE). */


boolean  InquireRSTDesc(RSTREE   rst,
                        char     *name,
                        int      *numbofdim,
                        int      *sizedirentry,
                        int      *sizedataentry,
                        int      *sizeinfo,
                        int      *maxdirfanout,
                        int      *maxdatafanout,
                        int      *pagesize,
                        int      *numbofdirpages,
                        int      *numbofdatapages,
                        int      pagesperlevel[],
                        int      *numbofrecords,
                        int      *height,
                        boolean  *unique);

         /* InquireRSTDesc provides some useful information about the
            R*-tree referenced by rst.
            
            name: see CreateRST.
            
            numbofdim:
            contains the number of dimensions of the R*-tree referenced by
            rst, i.e. the value the constant NumbOfDim was set to when it
            was created.
            
            sizedirentry, sizedataentry:
            contain the size (in bytes) of a directory and data entry
            respectively.
            
            sizeinfo:
            contains the size (in bytes) of an information part.
            
            maxdirfanout, maxdatafanout:
            contain the maximum possible number of entries a directory
            and data node respectively can store.
            
            pagesize: see CreateRST.
            
            numbofdirpages, numbofdatapages:
            total number of directory and data pages respectively in use.
            
            pagesperlevel:
            For each level i, beginning at the root, pagesperlevel[i]
            contains the number of pages in use.
            
            numbofrecords:
            total number of data records stored in the R*-tree.
            
            height:
            height of the tree, the lowest height is "1" (only the root
            exists).
            
            unique: see CreateRST. */


boolean  InsertRecord(RSTREE   rst,
                      typrect  rectangle,
                      typinfo  *info,
                      boolean  *inserted);

         /* InsertRecord inserts a new record in the R*-tree.
            If the unique flag is set TRUE (see CreateRST) a new record is
            not inserted if a record with the same rectangle is found. In
            this case inserted yields FALSE, but the return value is TRUE
            (if no error occurred).
            
            rectangle:
            is the rectangle part of the new record.
            
            info:
            is the information part of the new record. */


boolean  DeleteRecord(RSTREE   rst,
                      typrect  rectangle,
                      boolean  *deleted);

         /* DeleteRecord deletes the first record with the given rectangle
            it finds. It provides a fast deletion in trees where entries
            are unique and may be used in trees where entries are not
            unique, to delete iteratively all entries with the same
            rectangle as passed. */


boolean  ExistsRegion(RSTREE         rst,
                      typrect        queryrect1,
                      typrect        queryrect2,
                      DirQueryProc   DirQuery,
                      DataQueryProc  DataQuery,
                      boolean        *recordfound);

         /* ExistsRegion performs a RegionQuery on the R*-tree referenced
            by rst. It stops after the first record satisfying the query
            condition is found.
            See also RegionQuery.
            
            recordfound:
            is set to TRUE if a record satisfying the query condition
            exists, otherwise FALSE. */


boolean  RegionCount(RSTREE         rst,
                     typrect        queryrect1,
                     typrect        queryrect2,
                     DirQueryProc   DirQuery,
                     DataQueryProc  DatQuery,
                     int            *recordcount);

         /* RegionCount performs a RegionQuery on the R*-tree referenced
            by rst. It does not return records but only counts the number
            of records found.
            See also RegionQuery.
            
            recordcount:
            is set to the number of records satisfying the query
            condition. */


boolean  RegionQuery(RSTREE           rst,
                     typrect          queryrect1,
                     typrect          queryrect2,
                     DirQueryProc     DirQuery,
                     DataQueryProc    DatQuery,
                     QueryManageProc  Manage,
                     void             *pointer);

         /* RegionQuery performs a RegionQuery on the R*-tree referenced
            by rst. Up to two query rectangles may be passed by queryrect1
            and queryrect2. Two different procedures have two be provided
            (DirQuery, DataQuery) which perform the query in the directory
            and the data level respectively. A third procedure (Manage) must
            be provided to deal with the records successively found.
            A query is closed either if it did not find an additional
            record satisfying the query condition or if the finish flag is
            set by the procedure Manage.
            See also DirQueryProc, DataQueryProc, QueryManageProc.
            
            queryrect1, queryrect2:
            query rectangles to be compared with directory rectangles and
            data rectangles respectively.
            
            DirQuery, DataQuery:
            Procedure parameters passing comparison procedures of type
            boolean.
            
            Manage:
            Procedure parameter passing a management procedure.
            Manage is called each time a new data rectangle satisfying the
            query condition is found.
            Procedures of type QueryManageProc may provide the following
            facilities:
            Inspection of the data records rectangle and info part.
            Communication to another structure, pointed to by pointer.
            To modify the info part (the rectangle cannot be modified), and
            label the node to be written back.
            To finish the query.
            
            pointer:
            Arbitrary pointer passed through to the procedure Manage. */

boolean RegionQueryInfo(RSTREE R,
			Check includes,
			Check intersects,
			void * data,
			typrect rect,
			typdirinfo * info);

boolean  AllQuery(RSTREE           rst,
                  QueryManageProc  Manage,
                  void             *pointer);

         /* AllQuery performs a fast query which returns all records stored
            in the R*-tree referenced by rst.
            
            Manage:
            Procedure parameter passing a management procedure.
            Manage is called each time a new data rectangle satisfying the
            query condition is found.
            The type QueryManageProc provides the following functions:
            To inspect a data records rectangle and info part.
            To copy a record to the location pointer points to.
            To modify the info part (the rectangle cannot be modified),
            and label the node to be written back.
            To finish the query.
            
            pointer:
            Arbitrary pointer passed through to the procedure Manage.
            
            Since AllQuery is designed to be fast it does not support the
            complete counting facility.

            See also RegionQuery. */

boolean  Update(RSTREE           rst);

         /* Updates the directory nodes typdirinfo. */

boolean  JoinCountNv(RSTREE         rst1,
                     RSTREE         rst2,
                     typrect        R1queryrect1,
                     typrect        R1queryrect2,
                     typrect        R2queryrect1,
                     typrect        R2queryrect2,
                     DirQueryProc   Dir1Query,
                     DataQueryProc  Data1Query,
                     DirQueryProc   Dir2Query,
                     DataQueryProc  Data2Query,
                     DirQueryProc   DirJoin,
                     DataQueryProc  DataJoin,
                     int            *paircount);

         /* JoinCountNv performs a Join on the two R*-trees referenced by
            rst1 and rst2. It does not return record pairs but only counts
            the number of record pairs found.
            See also JoinNv.
            
            paircount:
            is set to the number of recordpairs satisfying the join
            condition. */


boolean  JoinNv(RSTREE          rst1,
                RSTREE          rst2,
                typrect         rst1queryrect1,
                typrect         rst1queryrect2,
                typrect         rst2queryrect1,
                typrect         rst2queryrect2,
                DirQueryProc    Dir1Query,
                DataQueryProc   Data1Query,
                DirQueryProc    Dir2Query,
                DataQueryProc   Data2Query,
                DirQueryProc    DirJoin,
                DataQueryProc   DataJoin,
                JoinManageProc  Manage,
                void            *pointer1,
                void            *pointer2);

         /* The functionality of the join can be considered as follows:
            On each of the two R*-trees to be joined a query is performed.
            On the resulting restricted sets of records of the two R*-trees
            the join is applied depending on the given join condition.
            
            Join performs a join on the two R*-trees referenced by rst1 and
            rst2. A join is closed either if it did not find an additional
            pair of records satisfying the join condition or if the finish
            flag is set by the procedure Manage.
            See also DirQueryProc, DataQueryProc, JoinManageProc.
            
            rst1queryrect1, rst1queryrect2:
            Used in connection with the pre-query on rst1;
            query rectangles to be compared with directory and data
            rectangles respectively of rst1.
            
            rst2queryrect1, rst2queryrect2:
            Used in connection with the pre-query on rst2;
            query rectangles to be compared with directory and data
            rectangles respectively of rst2.
            
            Dir1Query, Data1Query:
            Procedure parameters passing comparison procedures of type
            boolean.
            The two procedures have to perform the pre-query on rst1.
            See also DirQueryProc, DataQueryProc and RegionQuery.
            
            Dir2Query, Data2Query:
            Procedure parameters passing comparison procedures of type
            boolean.
            The two procedures have to perform the pre-query on rst2.
            See also DirQueryProc, DataQueryProc and RegionQuery.
            
            DirJoin, DataJoin:
            Procedure parameters passing comparison procedures of type
            boolean.
            DirJoin has to determine if two directory rectangles, one of
            rst1 the other of rst2 have to be joined,
            DataJoin has to determine if two data rectangles, one of
            rst1 the other of rst2 have to be joined.
            Here the types DirQueryProc and DataQueryProc are used as
            follows:
            DQP(RSTREE   rst1,
                typrect  rst1rect,
                typrect  rst2rect,
                typrect  unused);
            rst1      contains rst1.
            rst1rect  contains a directory or data rectangle respectively
                      of rst1.
            rst2rect  contains a directory or data rectangle respectively
                      of rst2.
            The last parameter is unused.
            
            Manage:
            Procedure parameter passing a management procedure.
            Manage is called each time a new pair of data rectangles
            satisfying the join condition is found.
            Procedures of type JoinManageProc may provide the following
            facilities:
            Inspection of the data records' rectangles and info parts.
            Communication to two other structures, pointed to by pointer1
            and pointer2.
            To finish the join.
            
            pointer1, pointer2:
            Arbitrary pointers passed through to the procedure Manage. */


/*************** ----- Performance Controll Routines ----- ***************/

/* ------ Counts-Switch: */

boolean  CountsOn0(RSTREE rst);
         /* put ON, set 0 */

boolean  CountsOn(RSTREE rst);
         /* put ON */

boolean  CountsOff(RSTREE rst);
         /* put OFF */
         
         /* the Counts-Switch applies to the variables set by the
            procedures GetCountRead, GetCountWrite.
            
            The procedure OpenRST initializes counting:
            the count variables are set to 0, the count switch is set to
            OFF. */


boolean  GetCountRead(RSTREE  rst,
                      int     *DirVisitCount,
                      int     *DataVisitCount,
                      int     *DirReadCount,
                      int     *DataReadCount);

         /* DirVisitCount is set to the number of directory nodes visited
            traversing the tree.
            DataVisitCount is set to the number of data nodes visited
            traversing the tree.
            DirReadCount directory nodes and DataReadCount data nodes had
            actually to be read from secondary memory.
            
            If the function returns FALSE these variables are set to 0.
            DirVisitCount, DataVisitCount is counted:
            for all query, join and update procedures.
            DirReadCount, DataReadCount is counted:
            whenever a read occurs. */


boolean  GetCountWrite(RSTREE  rst,
                       int     *DirModifyCount,
                       int     *DataModifyCount,
                       int     *DirWriteCount,
                       int     *DataWriteCount);

         /* DirModifyCount is set to the number of directory nodes
            modified.    
            DataModifyCount is set to the number of data nodes modified.
            DirWriteCount is set to the number of directory nodes written
            to secondary memory.
            DataWriteCount is set to the number of data nodes written to
            secondary memory.
            If the function returns FALSE these variables are set to 0.
            DirModifyCount, DataModifyCount is counted:
            for all update procedures.
            DirWriteCount, DataWriteCount is counted:
            whenever a write occurs. */


boolean  GetMemory(RSTREE  rst,
                   int     *numbofdirpages,
                   int     *numbofdatapages);

         /* numbofdirpages and numbofdatapages are set to the number of
            pages in use for the directory and the data level respectively.
            If the function returns FALSE these variables are set to 0. */


boolean  GetHeight(RSTREE  rst,
                   int     *height);

         /* height is set to the height of the tree, the lowest height
            is "1" (only the root exists).
            If the function returns FALSE these variables are set to 0. */


/*************************************************************************/

/* Layout of a directory and data node respectively (pidgin C)
   -----------------------------------------------------------
   
   directory node layout:
   
   struct {
          typrect  rectangle;
          int      subtree_pointer;
          } directory_entry;
   
   struct {
          int              n;
          directory_entry  entries[M];
          } directory_node;
   
   The maximum number of entries M varies with different page sizes. The
   minimum is M = 3.
   The minimum number of entries m is calculated as max(round(0.4 * M), 2).
   
   data node layout:
   
   struct {
          typrect  rectangle;
          typinfo  information_part;
          } data_entry;
   
   struct {
          int         n;
          data_entry  entries[M];
          } data_node;
   
   The maximum number of entries M varies with different page sizes. The
   minimum is M = 1.
   The minimum number of entries m is calculated as max(round(0.4 * M), 1).
*/

/*************************************************************************/

/* BUGS:

   The implementation does not provide packing, thus depending on the
   design of the machine and the compiler, and the choice of typatomkey and
   typinfo nodes may have gaps, i.e. the fanout may be smaller than you
   expect.
   Alignment problems, i.e. gaps between the entities stored in the nodes
   cause warning messages on stdout (with one restriction):   
   The implementation does not know the internal structure of typinfo (it
   only knows its size). Thus, if typinfo intrinsically contains gaps, this
   cannot be detected. Information about the actual values of important
   parameters may be obtained by calling InquireRSTDesc.
   
   R*-tree identifiers are only checked for NULL and non NULL. Though
   passing a non NULL invalid identifier is not detected as an error.
   
   To open an R*-tree twice may damage consistency. But the implementation
   does not detect this mistake.
   
   Early detection of memory limitations is not available. If accidentally
   the file system fills up during an update operation, the tree may be
   left in an inconsistent state.
   
   Although (aside from the lacks mentioned above) internal error detection
   is nearly exhaustive, all you get is a boolean return value.
   
   Passing the same R*-tree identifier twice to the JoinXX procedures is
   save. If a join is performed on one and the same R*-tree, a second
   R*-tree is opened virtually, but the performance control parameters are
   only available for one of them. Since the join internally does not work
   symmetrically this information generally is useless.
   
   If deletions are performed the files holding the directory and data
   pages will not shrink. Free pages are reoccupied by following
   insertions. A file reorganization algorithm is not implemented.
   
   The procedures of type QueryManageProc and JoinManageProc permit
   unprotected access to parts of the internal data structure of the
   R*-tree. This avoids copying but is unsafe of course.
   
   The implementation restricts the informational part of a data record to
   contain at least an integer sized contents.

*/


/*************************************************************************/
/*************** ----- For Private Test Purpose Only ----- ***************/

boolean  Find(RSTREE   rst,
              typrect  rectangle,
              boolean  *found,
              void     *buf,
              int      nbytes);


#endif /* !__RStarTree_h */
