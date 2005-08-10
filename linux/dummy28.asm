; Copyright 1995-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements creates some unused data items to get
; segments aligned on 32-byte boundaries in Linux.
;

	TITLE   setup

	.386

_GWDATA SEGMENT PAGE PUBLIC 'DATA'
dummy1		DD	0
dummy2		DD	0
dummy3		DD	0
dummy4		DD	0
dummy5		DD	0
dummy6		DD	0
dummy7		DD	0
_GWDATA ENDS

END
