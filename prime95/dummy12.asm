; Copyright 1995-1999 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements creates some unused data items to get
; the segment aligned on 32-byte boundary in Linux.
;

	TITLE   setup

	.386

_TEXT32A SEGMENT PAGE USE32 PUBLIC 'DATA'

dummy1		DD	0
dummy2		DD	0
dummy3		DD	0

_TEXT32A ENDS

END
