/*
 * Primenet communication routines for all operating systems
 * Uses sockets and HTTP
 */ 

/*
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//
// Copyright (c) 1997-2004 Just for Fun Software, Entropia Inc. and
//		 Peter Hunter.
// Used by Premission Only for Great Internet Mersenne Prime Search.
// All Rights Reserved.
//
//  MODULE:   primenet.c
//
//  PURPOSE:  Implements PrimeNet Version 4 API, as HTTP network client
//
//  AUTHOR:   Peter Hunter, on the basis of work by Scott Kurowski (v3 API)
//            Michiel van Loon, OS/2 adaptations 
//            Kurowski 5/1998, 4.0 API support for MPrime 16.x
//            Kurowski 9/1999, 4.0 API changes for MPrime 19.x
//	      Woltman 1/2002, Windows support and bug fixes
//          
//
//  ASSUMPTIONS: 1. less than 1k of data is sent or received per call
//               2. HTTP/1.0
//               3. PrimeNet Version 4 or later API on server and client
*/

/* Linux defines, adapted for OS/2, FreeBSD, and Windows */

#define PRIMENET_NO_ERROR 0

#ifdef __WATCOMC__
#include <types.h>
#endif
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>

#ifdef _WINDOWS_
#include <winsock.h>
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <errno.h>
#ifndef __IBMC__
#include <arpa/inet.h>
#include <unistd.h>
#endif
#endif

#if defined (__EMX__) && ! defined (AOUT)
#define htonl(x) (lswap(x))
#define ntohl(x) (lswap(x))
#define htons(x) (bswap(x))
#define ntohs(x) (bswap(x))
unsigned short bswap(unsigned short);
#endif

/* implement the missing inet_aton call */

#if defined (_WINDOWS_) || defined (__EMX__) || defined (__IBMC__)
int inet_aton (char *cp, struct in_addr *inp)
{
	u_long temp;
 
	temp = inet_addr(cp);
	if (temp == -1) return (0);
	inp->s_addr = temp;
	return (1);
}
#endif

/* Routine to get the error number after a failed sockets call */

int getLastSocketError (void)
{
#ifdef _WINDOWS_
	return (WSAGetLastError ());
#else
	return (errno);
#endif
}

/* simple password de/scrambler */

char SCRAMBLE_STRING[] = "/cgi-bin/pnHttp.exe";
char hx[] = "0123456789ABCDEF";

void scramble (char *s)
{
	char	out[100];
	char	*p = s, *z = out;
	unsigned int i, c = strlen (SCRAMBLE_STRING);

	for (i = 0; i < strlen (s); i++) {
		int b = (unsigned char) *p++ ^ SCRAMBLE_STRING[i % c];
		*z++ = hx[b >> 4];
		*z++ = hx[b % 16];
	}
	*z = 0;
	strcpy (s, out);
}

void unscramble (char *s)
{
	char	out[50];
	char	*q = s, *z = out;
	unsigned int i, c = strlen (SCRAMBLE_STRING);

	for (i = 0; i < strlen (s) >> 1; i++) {
		*z = (strchr (hx, *q++) - hx) * 16;
		*z += (strchr (hx, *q++) - hx);
		*z++ ^= SCRAMBLE_STRING[i % c];
	}
	*z = 0;
	strcpy (s, out);
}

/* base64 encode for basic proxy passwords */

static int encode[] = {
  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
  'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
  'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
  'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
  'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
  'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
  'w', 'x', 'y', 'z', '0', '1', '2', '3',
  '4', '5', '6', '7', '8', '9', '+', '/'
};

/* Base64-encode a null-terminated string. */

void encode64 (
	char	*data)
{
	char	*s, *end, *buf;
	unsigned int x, length;
	int	i, j;
	char	temp[400];

	length = strlen (data);
	if (length == 0) return;

	end = data + length - 3;

	buf = temp;
	for (s = data; s < end; ) {
		x = *s++ << 24;
		x |= *s++ << 16;
		x |= *s++ << 8;

		*buf++ = encode[x >> 26];
		x <<= 6;
		*buf++ = encode[x >> 26];
		x <<= 6;
		*buf++ = encode[x >> 26];
		x <<= 6;
		*buf++ = encode[x >> 26];
	}
	end += 3;

	x = 0;
	for (i = 0; s < end; i++)
		x |= *s++ << (24 - 8 * i);

	for (j = 0; j < 4; j++) {
		if (8 * i >= 6 * j) {
			*buf++ = encode [x >> 26];
			x <<= 6;
		} else {
			*buf++ = '=';
		}
	}

	*buf = 0;

	strcpy (data, temp);
}


/*///////////////////////////////////////////////////////////////////////////
//
// HTTP POST and GET procedure
//
///////////////////////////////////////////////////////////////////////////*/

#define INIFILENAME "primenet.ini"

char szSITE[] = "mersenne.org";		/* PrimeNet Server's home domain */
#define nHostPort 80			/* Internet PrimeNet port */
char szFILE[] = "/cgi-bin/pnHttp.exe";	/* interface application for GET */

/*
// pnHttpServer: POSTs or uses GET to send a formatted HTTP argument string
//               and downloads the server result page
*/

int pnHttpServer (char *pbuf, unsigned cbuf, char* postargs)
{
	char	szBuffer[1000];
	char	szURL[1024];			/* URL assembly buffer */
	char	buf[1024];
	struct in_addr defaddr;
	struct hostent *hp, def;
	struct sockaddr_in sn;
	int	s, res, debug, url_format;
	char	*alist[1];
	unsigned int count;
	char	szProxyHost[120], szUser[50], szPass[50], *con_host;
	unsigned short nProxyPort, con_port;

/* Get debug logging and URL format flags */

	debug = IniGetInt (INIFILENAME, "Debug", 0);
	url_format = IniGetInt (INIFILENAME, "UseFullURL", 2);
 
/* Get the host name of the optional proxy server.  If using a proxy */
/* server strip the optional http:// prefix. */

	IniGetString (INIFILENAME, "ProxyHost",
		      szProxyHost, sizeof (szProxyHost), NULL);
	if (*szProxyHost) {
		char	*colon;

		if ((szProxyHost[0] == 'H' || szProxyHost[0] == 'h') &&
		    (szProxyHost[1] == 'T' || szProxyHost[1] == 't') &&
		    (szProxyHost[2] == 'T' || szProxyHost[2] == 't') &&
		    (szProxyHost[3] == 'P' || szProxyHost[3] == 'p') &&
		    szProxyHost[4] == ':' && szProxyHost[5] == '/' &&
		    szProxyHost[6] == '/')
			strcpy (szProxyHost, szProxyHost + 7);

/* Get optional port number */

		if ((colon = strchr (szProxyHost, ':'))) {
			nProxyPort = (unsigned short) atoi (colon + 1);
			*colon = 0;
		} else
			nProxyPort = 8080;
		con_host = szProxyHost;
		con_port = nProxyPort;

/* Secure proxy - get username and password to negotiate access */

		IniGetString (INIFILENAME, "ProxyUser",
			      szUser, sizeof (szUser) - 1, NULL);
		IniGetString (INIFILENAME, "ProxyPass",
			      szPass, sizeof (szPass) - 1, NULL);

/* Scramble or unscramble the password as necessary */

		if (!IniGetInt (INIFILENAME, "ProxyMask", 0)) {
			scramble (szPass);
			IniWriteString (INIFILENAME, "ProxyPass", szPass);
			IniWriteInt (INIFILENAME, "ProxyMask", 1);
		}
		unscramble (szPass);
	}

/* No proxy server - use default site (mersenne.org) */

	else {
		con_host = szSITE;
		con_port = nHostPort;
	}

/* Output debug info */

redirect:
	if (debug) {
		sprintf (buf, "host = %s, port = %d\n", con_host, con_port);
		LogMsg (buf);
	}

/* Convert host name into an IP address */

	hp = gethostbyname (con_host);
	if (!hp) {
		char	szAltSiteAddr[30];
		if (con_host == szSITE) {
			IniGetString (INIFILENAME, "EntropiaIP",
				      szAltSiteAddr, 29, NULL);
			con_host = szAltSiteAddr;
		}
		if (!inet_aton (con_host, &defaddr)) {
			if (debug) {
				sprintf (buf, "PrimeNet: unknown host: %s\n",
					 con_host);
				LogMsg (buf);
			}
			return (PRIMENET_ERROR_CONNECT_FAILED);
		}
		alist[0] = (char *) &defaddr;
		def.h_name = szSITE;
		def.h_addr_list = alist;
		def.h_length = sizeof (struct in_addr);
		def.h_addrtype = AF_INET;
		def.h_aliases = 0;
		hp = &def;
	}

	if (debug) {
		sprintf (buf, "IP-addr = %s\n",
			 inet_ntoa (*(struct in_addr *)hp->h_addr));
		LogMsg (buf);
	}

	memset (&sn, 0, sizeof (sn));
	sn.sin_family = hp->h_addrtype;
	if (hp->h_length > (int) sizeof (sn.sin_addr)) {
		hp->h_length = sizeof (sn.sin_addr);
	}
	memcpy (&sn.sin_addr, hp->h_addr, hp->h_length);
	sn.sin_port = htons (con_port);

/* Create a socket and connect to server */

rel_url:
	if ((s = socket (hp->h_addrtype, SOCK_STREAM, 0)) < 0) {
		if (debug) {
			sprintf (buf, "Error in socket call: %d\n",
				 getLastSocketError ());
			LogMsg (buf);
		}
		return (PRIMENET_ERROR_CONNECT_FAILED);
	}

	if (connect (s, (struct sockaddr *) &sn, sizeof (sn)) < 0) {
		if (debug) {
			sprintf (buf, "Error in connect call: %d\n",
				 getLastSocketError ());
			LogMsg (buf);
		}
		closesocket (s);
		return (PRIMENET_ERROR_CONNECT_FAILED);
	}

/* GET method, data follows ? in URL */

	strcpy (szURL, "GET ");
	if (*szProxyHost || url_format) {
		strcat (szURL, "http://");
		strcat (szURL, szSITE);
		strcat (szURL, ":");
		sprintf (szBuffer, "%d", nHostPort);
		strcat (szURL, szBuffer);
	}
	strcat (szURL, szFILE);
	strcat (szURL, "?");
	strcat (szURL, postargs);
	strcat (szURL, " HTTP/1.0\r\n");

/* Append proxy authorization here */

	if (*szProxyHost && *szUser) {
		char	buf[200];
		strcat (szURL, "Proxy-Authorization: Basic ");
		sprintf (buf, "%s:%s", szUser, szPass);
		encode64 (buf);
		strcat (szURL, buf);
		strcat (szURL, "\r\n");
	}

	strcat (szURL, "\r\n");
	if (debug) LogMsg (szURL);

/* Send the URL request */

	res = send (s, szURL, strlen (szURL), 0);
	if (res < 0) {
		if (debug) {
			sprintf (buf, "Error in send call: %d\n",
				 getLastSocketError ());
			LogMsg (buf);
		}
		closesocket (s);
		return (PRIMENET_ERROR_SEND_FAILED);
	}

/* Now accumulate the response */

	*pbuf = 0; count = 1;
	while (count < cbuf) {
		res = recv (s, szBuffer, 999, 0);
		if (res < 0) {
			if (debug) {
				sprintf (buf, "Error in recv call: %d\n",
					 getLastSocketError ());
				LogMsg (buf);
			}
			closesocket (s);
			return (PRIMENET_ERROR_RECV_FAILED);
		}
		if (res == 0) break;
		szBuffer[res] = 0;
		if (debug) {
			sprintf (buf, "RECV: %s\n", szBuffer);
			LogMsg (buf);
		}
		if ((count + strlen(szBuffer)) >= cbuf)
			szBuffer[cbuf - count] = 0;
		strcat (pbuf, szBuffer);
		count += strlen (szBuffer);
	}

	closesocket (s);

/* pbuf + 9 is where the message code following HTTP/1.0 starts */

	if (count <= 10) res = -1;
	else res = atoi (pbuf + 9);

/* Some proxy servers can redirect us to another host.  These are */
/* the 300 series of error codes.  This code probably isn't right. */
/* We can improve it as people find problems. */

	if (res >= 300 && res <= 399) {
		char	*location, *colon;
		location = strstr (pbuf, "Location:");
		if (location != NULL) {
			if (debug) LogMsg ("Attempting a redirect.\n");
			location += 9;
			while (isspace (*location)) location++;
			strcpy (szProxyHost, location);

/* Parse the redirection address */

			if ((location[0] == 'H' || location[0] == 'h') &&
			    (location[1] == 'T' || location[1] == 't') &&
			    (location[2] == 'T' || location[2] == 't') &&
			    (location[3] == 'P' || location[3] == 'p') &&
			    location[4] == ':' && location[5] == '/' &&
			    location[6] == '/')
				strcpy (location, location + 7);
			con_host = location;

/* Get optional port number */

			if ((colon = strchr (location, ':')) != NULL) {
				con_port = (unsigned short) atoi (colon + 1);
				*colon = 0;
			} else
				con_port = 80;
			goto redirect;
		}
	}

/* Any return code other than 200 is an error.  We've had problems using */
/* both full URLs and relative URLs.  Thus, our default behavior is to try */
/* both before giving up. */

	if (res != 200) {
		if (url_format == 2 && *szProxyHost == 0) {
			if (debug) LogMsg ("Trying relative URL\n");
			url_format = 0;
			goto rel_url;
		}
		if (debug) LogMsg ("Return code is not 200\n");
		return (PRIMENET_ERROR_SERVER_UNSPEC);
	}

	return (PRIMENET_NO_ERROR);
}


/*///////////////////////////////////////////////////////////////////////////
//
// HTTP POST/GET argument formatting procedures
//
////////////////////////////////////////////////////////////////////////////*/

/* armor parameter control chars as hex codes for transport */

#define ARMOR_CHARS		"&+%\r\n"

void armor (char *s)
{
	char	buf[1024];
	char	*z = buf, *x = s;

/* & is token delimiter, '+' is space char */

	while (*x) {
		if (strchr (ARMOR_CHARS, *x)) {	
			*z++ = '%';	/* convert chars to %nn hex codes */
			*z++ = hx[(*x) / 16];
			*z++ = hx[(*x) % 16];
		} else if (*x == ' ')	/* convert spaces to '+' */
			*z++ = '+';
		else *z++ = *x;		/* copy normal character */
		x++;
	}
	*z = 0;
	strcpy (s, buf);
}

/* validate a text field against embedded control chars */
/* and optionally against spaces */

int validate (char *s, int fnospaces)
{ 
	/* illegal : &, %, + */

	while (*s) {
		if (*s == '&' || *s == '%' || *s == '+' ||
		    (fnospaces && *s == ' ')) return (1);
		s++;
	}
	return (0);
}


/* set a place-holding '.' character for a null string value */

void setifnull (char *s)
{
	if (*s) return;
	*s++ = '.';
	*s = 0;
}


/* validate a userID, userPW and computerID combination */

int validateID (void *pkt)
{
	struct primenetUserInfo *z = (struct primenetUserInfo*) pkt;

	if (validate (z->userID, 1)) return (1);
	if (validate (z->userPW, 1)) return (1);
	if (validate (z->computerID, 1)) return (1);

	setifnull (z->userID);
	setifnull (z->userPW);
	setifnull (z->computerID);

	return (0);
}

/*
// format_args: format a HTTP POST/GET argument string
//              from a PrimeNet 3.1 v4 packet
*/

int format_args (char* args, short operation, void* pkt)
{
	*args = 0;

	switch (operation) {
	case PRIMENET_PING_SERVER_INFO:		/* server info request */
		{
		struct primenetPingServerInfo *z = (struct primenetPingServerInfo*) pkt;
		sprintf (args,"ps&%d&.&.", z->u.serverInfo.versionNumber);
		break;
		}
	case PRIMENET_MAINTAIN_USER_INFO:	/* update/add user account */
		{
		struct primenetUserInfo *z = (struct primenetUserInfo*) pkt;
		if (validateID (z) ||
		    validate (z->oldUserID, 1) ||
		    validate (z->oldUserPW, 1))
			return (PRIMENET_ERROR_INVALID_PARAM);

		armor (z->userName);
		armor (z->userEmailAddr);

		setifnull (z->userName);
		setifnull (z->userEmailAddr);
		setifnull (z->oldUserID);
		setifnull (z->oldUserPW);

		sprintf (args,"uu&%d&%d&%d&%u&%u&%s&%s&%s&%s&%s&%s&%s&%d",
			 z->versionNumber, z->structSize, z->hashFunction,
			 z->hash, z->salt, z->userID, z->userPW,
			 z->computerID, z->userName, z->userEmailAddr,
			 z->oldUserID, z->oldUserPW, z->bUserOptions);
		break;
		}
	case PRIMENET_GET_ASSIGNMENT:		/* get assignment */
		{
		struct primenetGetAssignment *z = (struct primenetGetAssignment*) pkt;
		if (validateID (z)) return (PRIMENET_ERROR_INVALID_PARAM);
		sprintf (args,"ga&%d&%d&%d&%u&%u&%s&%s&%s&%d&%d&%04.1f",
			 z->versionNumber, z->structSize, z->hashFunction,
			 z->hash, z->salt, z->userID, z->userPW, z->computerID,
			 z->requestType, z->programType, z->how_far_factored);
		break;
		}
	case PRIMENET_COMPLETION_DATE:		/* update completion date */
		{
		struct primenetCompletionDate *z = (struct primenetCompletionDate*) pkt;
		if (validateID (z)) return (PRIMENET_ERROR_INVALID_PARAM);
		sprintf (args,"cd&%d&%d&%d&%u&%u&%s&%s&%s&%ld&%ld&%d&%d&%ld&%ld",
			 z->versionNumber, z->structSize, z->hashFunction,
			 z->hash, z->salt, z->userID, z->userPW, z->computerID,
			 z->exponent, z->days, z->requestType, z->programType,
			 z->iteration, z->nextMsg);
		break;
		}
	case PRIMENET_RESULT_MESSAGE:		/* log a text message */
		{
		struct primenetResultMessage *z = (struct primenetResultMessage*) pkt;
		if (validateID (z)) return (PRIMENET_ERROR_INVALID_PARAM);

		armor (z->message);

		sprintf (args,"rm&%d&%d&%d&%u&%u&%s&%s&%s&%s",
			 z->versionNumber, z->structSize, z->hashFunction,
			 z->hash, z->salt, z->userID, z->userPW, z->computerID,
			 z->message);
		break;
		}
	case PRIMENET_ASSIGNMENT_RESULT:	/* assignment result */
		{
		struct primenetAssignmentResult *z = (struct primenetAssignmentResult*) pkt;
		if (validateID (z)) return (PRIMENET_ERROR_INVALID_PARAM);
		if (z->resultType == PRIMENET_RESULT_NOFACTOR)
			sprintf (args, "ar&%d&%d&%d&%u&%u&%s&%s&%s&%ld&%d&%04.1f",
				 z->versionNumber, z->structSize,
				 z->hashFunction, z->hash, z->salt, z->userID,
				 z->userPW, z->computerID, z->exponent,
				 z->resultType,
				 z->resultInfo.how_far_factored);
		else if (z->resultType == PRIMENET_RESULT_FACTOR ||
			 z->resultType == PRIMENET_RESULT_TEST ||
			 z->resultType == PRIMENET_RESULT_PRIME)
			sprintf (args, "ar&%d&%d&%d&%u&%u&%s&%s&%s&%ld&%d&%s",
				 z->versionNumber, z->structSize,
				 z->hashFunction, z->hash, z->salt, z->userID,
				 z->userPW, z->computerID, z->exponent,
				 z->resultType, z->resultInfo.factor);
				/* all return szString[] from same offset */
		else if (z->resultType == PRIMENET_RESULT_UNRESERVE)
			sprintf (args, "ar&%d&%d&%d&%u&%u&%s&%s&%s&%ld&%d",
				 z->versionNumber, z->structSize,
				 z->hashFunction, z->hash, z->salt, z->userID,
				 z->userPW, z->computerID, z->exponent,
				 z->resultType);
		else
			return (PRIMENET_ERROR_INVALID_PARAM);
		break;
		}
	case PRIMENET_SET_COMPUTER_INFO: /* set client's machine information */
		{
		struct primenetComputerInfo *z = (struct primenetComputerInfo*) pkt;
		if (validateID (z)) return (PRIMENET_ERROR_INVALID_PARAM);
		sprintf (args,"mi&%d&%d&%d&%u&%u&%s&%s&%s&%d&%d&%d",
			 z->versionNumber, z->structSize, z->hashFunction,
			 z->hash, z->salt, z->userID, z->userPW, z->computerID,
			 z->cpu_type, z->speed, z->hours_per_day);
		break;
		}
	default:
		return (PRIMENET_ERROR_INVALID_FIELD);	/* default error */
	}
	return (0);
}


/*////////////////////////////////////////////////////////////////////////////
//
// HTTP downloaded response page parsing procedures
//
/////////////////////////////////////////////////////////////////////////////*/

/* skip over the token name and point to the data string */

char* skip_to_pnResult (char *s)
{
	char	*p;
	p = strstr (s, "pnResult=");
	if (p == NULL) return (s + strlen (s));
	return (p + 9);
}

char* skip_token (char *s)
{
	while (*s && *s != '=') s++;
	if (*s == '=') s++;
	return (s);
}


/* copy the data string up to the next '\r' delimiter character */

char* copy_value (char *buf, char *s)
{
	while (*s && *s != '\r') *buf++ = *s++;
	if (*s == '\r') s++;
	*buf = 0;
	return (s);
}


/*
// parse_page: reads the server response page tokens and values
//             and converts these back into a PrimeNet 3.1 v4 packet
*/

int parse_page (char *buf, short operation, void *pkt)
{
	char	*s, *stop;
	char	lbuf[24];
	int	res;

/* get result code, which is always first */

	s = skip_to_pnResult (buf);
	if (*s == 0 || *s < '0' || *s > '9')
		return (PRIMENET_ERROR_HTTP_BAD_PAGE);

	res = strtol (s, &s, 10);

/* quit now if garbage (for some reason) follows result code */

	if (*s != '\r') return (res);

	switch (operation) {
		/* update completion, results message, assignment result */
		/* each have no data other than a result code */
	case PRIMENET_COMPLETION_DATE:
	case PRIMENET_ASSIGNMENT_RESULT:
	case PRIMENET_RESULT_MESSAGE:
	case PRIMENET_SET_COMPUTER_INFO:
		break;
	case PRIMENET_PING_SERVER_INFO:
		/* ping returns server info */
		{
		struct primenetPingServerInfo *z = (struct primenetPingServerInfo*) pkt;

		s = copy_value (lbuf, skip_token (s));
		z->u.serverInfo.versionNumber = (short) strtol (lbuf, &stop, 10);

/* v19+ message info ping packet */

		if (z->u.serverInfo.versionNumber >= PRIMENET_PING_INFO_MESSAGE) {
			s = copy_value (lbuf, skip_token (s));
			z->u.messageInfo.msgID = (short) strtol (lbuf, &stop, 10);
			s = copy_value (lbuf, skip_token (s));
			z->u.messageInfo.targetMap = strtoul (lbuf, &stop, 10);
			s = copy_value (z->u.messageInfo.szMsgText, skip_token (s));
		}

/* normal server info ping packet */

		else {
			s = copy_value (z->u.serverInfo.buildID, skip_token (s));
			s = copy_value (z->u.serverInfo.primenetServerName, skip_token (s));
			s = copy_value (z->u.serverInfo.adminEmailAddr, skip_token (s));
		}
		break;
		}
	case PRIMENET_MAINTAIN_USER_INFO:
		{
		/* update user info can return userID, userPW, */
		/* computerID (usually what's sent) */

		struct primenetUserInfo *z = (struct primenetUserInfo*) pkt;

		s = copy_value (z->userID, skip_token (s));
		s = copy_value (z->userPW, skip_token (s));
		s = copy_value (z->computerID, skip_token (s));

		if (z->userID[0] == '.') z->userID[0] = 0;
		if (z->userPW[0] == '.') z->userPW[0] = 0;
		if (z->computerID[0] == '.') z->computerID[0] = 0;

		break;
		}
	case PRIMENET_GET_ASSIGNMENT:
		{
		/* get assignment returns requestType, exponent, */
		/* how_far_factored */
	
		struct primenetGetAssignment *z = (struct primenetGetAssignment*) pkt;

		s = copy_value (lbuf, skip_token (s));
		z->requestType = (short) strtol (lbuf, &stop, 10);
		s = copy_value (lbuf, skip_token (s));
		z->exponent = strtol (lbuf, &stop, 10);
		s = copy_value (lbuf, skip_token (s));
		z->how_far_factored = strtod (lbuf, &stop);

		break;
		}
	default:
		return (PRIMENET_ERROR_HTTP_BAD_PAGE);
	}
	return (res);
}


/*
// Primenet: main interface to Prime95.exe, accepts
//           and returns PrimeNet 3.0 packets
*/

int PRIMENET (short operation, void *pkt)
{
	int	status;
	char args[1024];		/* formatted arguments buffer */
	char pbuf[4096];		/* return page buffer */

/* Assemble GET/POST arguments */

	status = format_args (args, operation, pkt);
	if (status != PRIMENET_NO_ERROR) return (status);

/* Send arguments, read back resulting page */

	status = pnHttpServer (pbuf, sizeof (pbuf), args);
	if (status != PRIMENET_NO_ERROR) return (status);

/* Extract results from returned page into packet */

	return (parse_page (pbuf, operation, pkt));
}
