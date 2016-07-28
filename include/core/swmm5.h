//-----------------------------------------------------------------------------
//   swmm5.h
//
//   Project: EPA SWMM5
//   Version: 5.1
//   Date:    03/24/14  (Build 5.1.001)
//   Author:  L. Rossman
//
//   Prototypes for SWMM5 functions exported to swmm5.dll.
//
//-----------------------------------------------------------------------------
#ifndef SWMM5_H
#define SWMM5_H


// --- define WINDOWS

#undef WINDOWS
#ifdef _WIN32
  #define WINDOWS
#endif
#ifdef __WIN32__
  #define WINDOWS
#endif

// --- define DLLEXPORT

#ifdef WINDOWS
#define DLLEXPORT __declspec(dllexport) __stdcall
#define SDLLEXPORT __declspec(dllexport) 
#define STDCALL __stdcall
#else
#define DLLEXPORT 
#define SDLLEXPORT
#define STDCALL
#pragma GCC visibility push(default)
#endif

// --- use "C" linkage for C++ programs

#ifdef __cplusplus
extern "C" { 
#endif 

typedef struct Project Project;
typedef double DateTime;
typedef struct TNode TNode;
typedef struct TLink TLink;
typedef struct TSubcatch TSubcatch;

int  DLLEXPORT   swmm_run(char* f1, char* f2, char* f3); 
SDLLEXPORT Project* STDCALL swmm_open(char* f1, char* f2, char* f3);
int  DLLEXPORT   swmm_start(Project* project, int saveFlag);
int  DLLEXPORT   swmm_step(Project* project, double* elapsedTime);
int  DLLEXPORT   swmm_end(Project* project);
int  DLLEXPORT   swmm_report(Project* project);
int  DLLEXPORT   swmm_getMassBalErr(Project* project, float* runoffErr, float* flowErr,float* qualErr);
int  DLLEXPORT   swmm_close(Project* project);
int  DLLEXPORT   swmm_getVersion();


//additional
int  DLLEXPORT   swmm_getErrorCode(Project* project);
double DLLEXPORT swmm_getDateTime(Project* project, char* beginorend);
void DLLEXPORT  datetime_decodeDateTime(DateTime dateTime, int* y, int* m, int* d, int* h, int* mm, int* s);
SDLLEXPORT char * getErrorMsg(int errorCode);
int DLLEXPORT getObjectTypeCount(Project* project, int type);

//TNode
SDLLEXPORT TNode* STDCALL getNode(Project* project, int index);
SDLLEXPORT TNode* STDCALL getNodeById(Project* project, char* id);
void DLLEXPORT setNode(Project* project, char* nodeId, char* propertyName, double value);

//TLink
SDLLEXPORT TLink*  STDCALL getLink(Project* project, int index);
SDLLEXPORT TLink* STDCALL getLinkById(Project* project, char* id);
void DLLEXPORT setLink(Project* project, char* linkId, char* propertyName, double value);

//TSubcatch
SDLLEXPORT TSubcatch*  STDCALL getSubcatch(Project* project, int index);
SDLLEXPORT TSubcatch* STDCALL getSubcatchById(Project* project, char* id);
void DLLEXPORT setSubcatch(Project* project, char* subCatchId, char* propertyName, double value);

#ifdef __cplusplus 
}   // matches the linkage specification from above */ 
#endif

#ifdef WINDOWS
#else
#pragma GCC visibility pop
#endif

#endif
