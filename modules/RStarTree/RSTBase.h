/* -----  RSTBase.h  ----- */
#ifndef __RSTBase_h
#define __RSTBase_h


/*** Begin    -----  naming rules for structured types  -----   **********/
/*
Non vector types are always called "typ<..>",
  Exceptions:
    Side, rstree.
References to non vector types are always called "ref<..>".
  Exceptions:

    RSTREE.
Vector types are always called <..>Array.
  Exceptions:
    RSTName, typpoint, typrect, typfixblock,
  and the internal arrays of the R*-tree as
    typpath, typentrynumbs, typpagenrs, typflagarray, typpagecountarr.
*/
/***   End    -----  naming rules for structured types  -----   **********/

/* constants */

#define SIZEfixblock 512
#define maxnum (((SIZEfixblock-2*sizeof(int)) / sizeof(int))-1)
#define maxentries MaxNumbOfEntriesPerPage
#define EntryRange (maxentries+sizeof(int))

#define chainlen 42     /* 3^ 40  =  12,157,665,459,056,928,801 */

#define datasuffix ".Data"
#define dirPDsuffix ".DirPD"
#define dataPDsuffix ".DataPD"

#define STDMODE 0644

#define paramblocknumb 0
#define firstPDblocknumb (paramblocknumb+1)
#define rootblocknumb 0


/* types */

typedef typinterval  *refinterval;

typedef char RSTName[MaxNameLength];
typedef double typpoint[NumbOfDim];

typedef double    ValueArray[EntryRange];
typedef int       IndexArray[EntryRange];
typedef byte      ByteArray[EntryRange];
typedef typrect   RectArray[EntryRange];
typedef typrect   NbrsArray[EntryRange];
typedef typpoint  PointArray[EntryRange];

typedef enum {
             low,
             high
             } Side;

typedef struct {
               typrect    rect;
               int        ptrtosub;
               PADDING_32_BITS;
               typdirinfo info;
               } typDIRent, *refDIRent;		/* inner entry */
typedef struct {
               typrect  rect;
               typinfo  info;
               } typDATAent, *refDATAent;	/* data entry */
typedef union {
              typDIRent   DIR;
              typDATAent  DATA;
              } typentry, *refentry;		/* universal entry */
	/* If typentry v; Then                  */
	/*   every access to rect by v.DIR.rect */
	/* End                                  */
	/* hope the compiler-designers agree.   */


typedef typDIRent   typDIRentries[EntryRange];
typedef typDATAent  typDATAentries[EntryRange];


typedef struct {
               int            nofentries;
               PADDING_32_BITS;
               typDIRentries  entries;		/* 0 .. nofentries-1 !! */
               } typDIRnode, *refDIRnode;	/* inner node */
typedef struct {
               int             nofentries;
               typDATAentries  entries;		/* 0 .. nofentries-1 !! */
               } typDATAnode, *refDATAnode;	/* data node */
typedef union {
              typDIRnode   DIR;
              typDATAnode  DATA;
              } typnode, *refnode;		/* universal node */
	/* If typnode v; Then                               */
	/*   every access to nofentries by v.DIR.nofentries */
	/* End                                              */
	/* hope the compiler-designers agree.               */


typedef refnode  typpath[chainlen+1];		/* 1 .. height !! */
typedef int      typentrynumbs[chainlen+1];	/* 1 .. height !! */
typedef int      typpagenrs[chainlen+1];	/* 1 .. height !! */
typedef boolean  typflagarray[chainlen+1];	/* 1 .. height !! */

typedef byte  typfixblock[SIZEfixblock];

typedef struct {
               int  childnr;
               int  nofnumbers;
               int  number[maxnum+1];		/* 0 + [1 .. maxnum] !! */
               } typpagedir, *refpagedir;

typedef union {
              typpagedir   _;			/* call the intrinsics "_" */
              typfixblock  fixblock;		/* adjust to SIZEfixblock   */
              } typPDblock, *refPDblock;

typedef int  typpagecountarr[chainlen+1];	/* 1 .. height !! */

typedef struct {
               boolean         unique;
               int             height;
               int             SIZE_DIRnofentries, SIZE_DATAnofentries;
               int             SIZEinfo;
               int             direntrylen, dataentrylen;
               int             nbrsexam;
               int             reinstpercent;
               int             minfillpercent;
               int             dirreinsertqty, datareinsertqty;
               int             pagelen;
               int             dirM, dirMwords, dirm;
               int             dataM, dataMwords, datam;
               int             maxdim;
               int             dirpagecount, datapagecount, recordcount;
               typpagecountarr pagecountarr;
               } typparameters, *refparameters;

typedef union {
              typparameters  _;			/* call the intrinsics "_" */
              typfixblock    fixblock;		/* adjust to SIZEfixblock   */
              } typparamblock, *refparamblock;

typedef struct { 
               boolean  countflag;
               int      dirvisitcount, datavisitcount;
               int      dirreadcount, datareadcount;
               int      dirmodifycount, datamodifycount;
               int      dirwritecount, datawritecount;
               } typcount, *refcount;

typedef struct {
               File f;
               int bl;
               } typfiledesc;

typedef struct {
               int            dirnodelen, datanodelen;
               typpath        N, NInst, NDel;
               typentrynumbs  E, EInst;
               typpagenrs     P;
               typflagarray   Nmodified, ReInsert;
               typPDblock     dirpagedir, datapagedir;
               typparamblock  parameters;
               typcount       count;
               boolean        RSTDone;
               PADDING_32_BITS;
               refnode        helpdirnode, helpdatanode, Ntosplit, Nsibling;
               typfiledesc    dir, data, dirPD, dataPD;
               RSTName        dirname;
               boolean        readonly;
               } rstree;


/****************************************************************/
/* begin **** --- types to check for alignment problems --- *****/

typedef
  typDIRent typ2DIRentries[2];
typedef
  typDATAent typ2DATAentries[2];
typedef
  struct {
           int nofentries;
           PADDING_32_BITS;
           typ2DIRentries entries;
         } typDIRnodeOf2;	/* inner node containing 2 entries */
typedef
  struct {
           int nofentries;
           typ2DATAentries entries;
         } typDATAnodeOf2;	/* data node containing 2 entries */

typedef
  typDIRent typ3DIRentries[3];
typedef
  typDATAent typ3DATAentries[3];
typedef
  struct {
           int nofentries;
           PADDING_32_BITS;
           typ3DIRentries entries;
         } typDIRnodeOf3;	/* inner node containing 3 entries */
typedef
  struct {
           int nofentries;
           typ3DATAentries entries;
         } typDATAnodeOf3;	/* data node containing 3 entries */

/*   end **** --- types to check for alignment problems --- *****/
/****************************************************************/


#endif /* !__RSTBase_h */
