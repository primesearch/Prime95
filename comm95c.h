/* Common definitions in Prime95, Saver95, and NTPrime */

extern HMODULE HLIB;			// Handle of networking DLL
extern int (__stdcall *PRIMENET)(short, void *);

/* No OS specific tasks to execute */

#define doMiscTasks()

