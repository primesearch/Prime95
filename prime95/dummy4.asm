; Copyright 1995-2002 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements creates some unused data items to get
; segments aligned on 32-byte boundaries in Linux.
;

	TITLE   setup

	.386

_DATA SEGMENT PAGE USE32 PUBLIC 'DATA'
dummy1		DD	0
_DATA ENDS

END
