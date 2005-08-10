/*----------------------------------------------------------------------
| common.h
|
| This file contains handy #defines that I use in all my projects
| 
|  Copyright 2005 Just For Fun Software, Inc.
|  All Rights Reserved.
+---------------------------------------------------------------------*/

#ifndef _COMMON_H
#define _COMMON_H

#ifndef TRUE
#define TRUE	1
#endif
#ifndef FALSE
#define FALSE	0
#endif

/* Define the ASSERT macro I use while debugging */

#ifdef GDEBUG
#include <assert.h>
#define ASSERTG		assert
#else
#define ASSERTG(a)
#endif

#endif
