/* #define _DEBUG */

#include "prime.h"

/* Linux defines */
/* Adapted for OS/2 and FreeBSD */

#define RPC_S_OK 0
typedef int error_status_t;

/*
 * Primenet routines for Linux and OS/2
 */ 

/*
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//
// Copyright (c) 1997-1999 Entropia.com, Inc. and Peter Hunter.
// Used by Premission Only for Great Internet Mersenne Prime Search.
// All Rights Reserved.
//
//  MODULE:   primenet.c (based on http_primenet.cpp)
//
//  PURPOSE:  Implements PrimeNet Version 4 API, as HTTP network client
//
//  AUTHOR:   Peter Hunter, on the basis of work by Scott Kurowski (v3 API)
//            Michiel van Loon, OS/2 adaptations 
//            Kurowski 5/1998, 4.0 API support for MPrime 16.x
//            Kurowski 9/1999, 4.0 API changes for MPrime 19.x
//          
//
//  ASSUMPTIONS: 1. less than 1k of data is sent or received per call
//               2. HTTP/1.0
//               3. PrimeNet Version 4 or later API on server and in calling client
*/

#if defined (__EMX__) && ! defined (AOUT)
#define htonl(x) (lswap(x))
#define ntohl(x) (lswap(x))
#define htons(x) (bswap(x))
#define ntohs(x) (bswap(x))
unsigned short bswap(unsigned short);
#endif

#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#ifndef __IBMC__
#include <arpa/inet.h>
#endif
#include <netdb.h>
#include <string.h>
#ifndef __IBMC__
#include <unistd.h>
#endif
#include <ctype.h>
#include <stdio.h>

#ifdef __EMX__
#define PRIMENET_TMP "primenet.tmp"
#endif

#if defined (__EMX__) || defined (__IBMC__)

int inet_aton(char *cp, struct in_addr *inp)
{
/* implement the missing inet_aton call */
  u_long temp;
 
  temp = inet_addr(cp);
  if (temp == -1) 
    return -1;
  else {
    inp->s_addr = temp;
    return 0;
  } 
}
#endif

/*///////////////////////////////////////////////////////////////////////////
//
// HTTP POST and GET procedure
//
///////////////////////////////////////////////////////////////////////////*/

#ifdef __linux__
#define INET_APP_NAME	"PrimeNet4Linux"	/* App name seen by server */
#endif
#ifdef __FreeBSD__
#define INET_APP_NAME	"PrimeNet4FreeBSD"	/* App name seen by server */
#endif
#ifdef __EMX__
#define INET_APP_NAME	"PrimeNet4OS2"		/* App name seen by server */
#endif

#define INIFILENAME "primenet.ini"

#define _USE_GET_METHOD		/* I am choosing GET since it skips the */
				/* intermediate perl script on the server */

char szSITE[] = "entropia.com";		/* PrimeNet Server's home domain */
int nHostPort = 80;			/* Internet PrimeNet port */

char szFILE[] = "/cgi-bin/pnHttp.exe";	/* interface application for GET */
char szURL[1024];			/* URL assembly buffer */

/*
// pnHttpServer: POSTs or uses GET to send a formatted HTTP argument string
//               and downloads the server result page
*/

char hx[] = "0123456789ABCDEF";

int pnHttpServer(char *pbuf, unsigned cbuf, char* postargs)
{
	char szBuffer[1000];
	struct in_addr defaddr;
	struct hostent *hp, def;
	struct sockaddr_in sn;
	int s;
#ifdef __EMX__
	int bytes_recv;
#endif
	FILE *fp;
	char *alist[1];
	unsigned int count;
	char szProxyHost[120], *con_host;
	char szAltSiteAddr[30] = "";
	int nProxyPort, con_port;
	char *colon;
	
	IniGetString(INIFILENAME, "ProxyHost", szProxyHost, 120, NULL);
	if (*szProxyHost) {
		if ((colon = strchr(szProxyHost, ':'))) {
			nProxyPort = atoi(colon + 1);
			*colon = 0;
		} else
			nProxyPort = 8080;
		con_host = szProxyHost;
		con_port = nProxyPort;
	} else {
		con_host = szSITE;
		con_port = nHostPort;
	}

#ifdef _DEBUG
	fprintf(stderr,"con_host = %s\n",con_host);
	fprintf(stderr,"con_port = %d\n",con_port);
#endif

	hp = gethostbyname(con_host);
	if (!hp) {
		if (con_host == szSITE) {
			IniGetString (INIFILENAME, "EntropiaIP", szAltSiteAddr, 29, NULL);
			con_host = szAltSiteAddr;
		}
		if (!inet_aton(con_host, &defaddr)) {
#ifdef _DEBUG
			fprintf(stderr, "PrimeNet: unknown host: %s\n", con_host);
#endif
			return PRIMENET_ERROR_CONNECT_FAILED;
		}
		alist[0] = (char *)&defaddr;
		def.h_name = szSITE;
		def.h_addr_list = alist;
		def.h_length = sizeof(struct in_addr);
		def.h_addrtype = AF_INET;
		def.h_aliases = 0;
		hp = &def;
	}
#ifdef _DEBUG
	fprintf(stderr,"IP-addr = %s\n",inet_ntoa(*(struct in_addr *)hp->h_addr));
#endif
	memset(&sn, 0, sizeof(sn));
	sn.sin_family = hp->h_addrtype;
	if (hp->h_length > (int)sizeof(sn.sin_addr)) {
		hp->h_length = sizeof(sn.sin_addr);
	}
	memcpy(&sn.sin_addr, hp->h_addr, hp->h_length);
	sn.sin_port = htons(con_port);
	if ((s = socket(hp->h_addrtype, SOCK_STREAM, 0)) < 0) {
#ifdef _DEBUG
	        perror("PrimeNet: socket");
#endif
		return PRIMENET_ERROR_CONNECT_FAILED;
	}

	if (connect(s, (struct sockaddr *)&sn, sizeof(sn)) < 0) {
#ifdef _DEBUG
		perror("PrimeNet: connect");
#endif
		close(s);
		return PRIMENET_ERROR_CONNECT_FAILED;
	}

	/* get method, data follows ? in URL */

	strcpy(szURL, "GET ");
	if (*szProxyHost) {
		strcat(szURL, "http://");
		strcat(szURL, szSITE);
		strcat(szURL, ":");
		sprintf(szBuffer, "%d", nHostPort);
		strcat(szURL, szBuffer);
	}
	strcat(szURL, szFILE);
	strcat(szURL, "?");
	strcat(szURL, postargs);
	strcat(szURL, " HTTP/1.0\r\n\r\n");
#ifdef _DEBUG
	puts(szURL);
#endif
#ifdef __EMX__
	send(s,szURL, strlen(szURL),0);
#else
	write(s, szURL, strlen(szURL));
#endif
	*szBuffer = 0;
#ifdef __EMX__
	fp = fopen(PRIMENET_TMP,"w+");
	while ((bytes_recv=recv(s,szBuffer,999,0)) >0) {
	  szBuffer[bytes_recv] = 0;
	  fwrite(szBuffer,1,bytes_recv,fp);
	}
	rewind(fp);
#else
	fp = fdopen(s, "r");
#endif	
        fgets(szBuffer, 999, fp);
#ifdef _DEBUG
	puts("First buffer");
	puts(szBuffer);
#endif
	
	/* szBuffer + 9 is where the message following HTTP/1.0 starts */
	/* TODO: more subtle checking here */

	if ((strlen(szBuffer)<=9) || (atoi(szBuffer + 9)!=200)) {
#ifdef _DEBUG
		fprintf(stderr, "Error from server: %s\n", szBuffer);
		fputs(szBuffer, stderr);
#endif
		fclose(fp);
		return PRIMENET_ERROR_SERVER_UNSPEC;
	}
#ifdef _DEBUG
	puts("Get succeeded!");
#endif

	count = 1; *pbuf = 0;
	while ((count < cbuf) && !feof(fp)) {
		*szBuffer = 0;
		fgets(szBuffer, 999, fp);
		if ((count + strlen(szBuffer)) >= cbuf)
			szBuffer[cbuf - count] = 0;
		strcat(pbuf, szBuffer);
		count += strlen(szBuffer);
#ifdef _DEBUG
		printf("%s",szBuffer);
#endif
	}
	fclose(fp);
#ifdef __EMX__
	unlink(PRIMENET_TMP);
#endif
	return RPC_S_OK;
}


/*///////////////////////////////////////////////////////////////////////////
//
// HTTP POST/GET argument formatting procedures
//
////////////////////////////////////////////////////////////////////////////*/

/* armor parameter control chars as hex codes for transport */

#define ARMOR_CHARS		"&+%\r\n"

void armor( char* s )
{
	char buf[1024];
	char *z = buf, *x = s;

	while ( *x )
	{
		if ( strchr( ARMOR_CHARS, *x ) )	/* & is token delimiter, '+' is space char */
		{
			*z++ = '%';			/* convert illegal chars to %nn hex codes */
			*z++ = hx[(*x) / 16];
			*z++ = hx[(*x) % 16];
		}
		else if ( *x == ' ' )			/* convert spaces to '+' for argument transport */
			*z++ = '+';
		else *z++ = *x;				/* copy normal character */
		x++;
	}
	*z = '\0';
	strcpy( s, buf );
}

/* validate a text field against embedded control chars, and optionally against spaces */

int validate( char* s, int fnospaces )
{ 
	/* illegal : &, %, + */

	while ( *s )
	{
		if (*s == '&' || *s == '%' || *s == '+' || (fnospaces && *s == ' ')) return 1;
		s++;
	}
	return 0;
}


/* set a place-holding '.' character for a null string value */

void setifnull( char* s )
{
	if ( *s ) return;
	*s++ = '.';
	*s = 0;
}


/* validate a userID, userPW and computerID combination */

int validateID( void* pkt )
{
	struct primenetUserInfo *z = (struct primenetUserInfo*) pkt;

	if ( validate( z->userID, 1 ) ) return 1;
	if ( validate( z->userPW, 1 ) ) return 1;
	if ( validate( z->computerID, 1 ) ) return 1;
	
	setifnull( z->userID );
	setifnull( z->userPW );
	setifnull( z->computerID );

	return 0;
}

/*
// format_args: format a HTTP POST/GET argument string
//              from a PrimeNet 3.1 v4 packet
*/

int format_args( char* args, short operation, void* pkt )
{
	*args = 0;

	switch ( operation )
	{
	case PRIMENET_PING_SERVER_INFO:		/* server info request */
		{
			struct primenetPingServerInfo *z = (struct primenetPingServerInfo*) pkt;

			sprintf(args,"ps&%d&.&.", z->u.serverInfo.versionNumber);
			break;
		}
	case PRIMENET_MAINTAIN_USER_INFO:	/* update/add user account */
		{
			struct primenetUserInfo *z = (struct primenetUserInfo*) pkt;

			if ( validateID( z ) ||
			     validate( z->oldUserID, 1 ) ||
				 validate( z->oldUserPW, 1 ) ||
				 validate( z->userEmailAddr, 1 ) )
				 return PRIMENET_ERROR_INVALID_PARAM;

			armor( z->userName );
			armor( z->userEmailAddr );

			setifnull( z->userName );
			setifnull( z->userEmailAddr );
			setifnull( z->oldUserID );
			setifnull( z->oldUserPW );

			sprintf(args,"uu&%d&%d&%d&%u&%u&%s&%s&%s&%s&%s&%s&%s&%d",
							z->versionNumber,
							z->structSize,
							z->hashFunction,
							z->hash,
							z->salt,
							z->userID,
							z->userPW,
							z->computerID,
							z->userName,
							z->userEmailAddr,
							z->oldUserID,
							z->oldUserPW,
							z->bUserOptions);
			break;
		}
	case PRIMENET_GET_ASSIGNMENT:		/* get assignment */
		{
			struct primenetGetAssignment *z = (struct primenetGetAssignment*) pkt;

			if ( validateID( z ) ) return PRIMENET_ERROR_INVALID_PARAM;

			sprintf(args,"ga&%d&%d&%d&%u&%u&%s&%s&%s&%d&%d&%04.1f",
							z->versionNumber,
							z->structSize,
							z->hashFunction,
							z->hash,
							z->salt,
							z->userID,
							z->userPW,
							z->computerID,
							z->requestType,
							z->programType,
							z->how_far_factored);
			break;
		}
	case PRIMENET_COMPLETION_DATE:		/* update completion date */
		{
			struct primenetCompletionDate *z = (struct primenetCompletionDate*) pkt;

			if ( validateID( z ) ) return PRIMENET_ERROR_INVALID_PARAM;

			sprintf(args,"cd&%d&%d&%d&%u&%u&%s&%s&%s&%ld&%ld&%d&%d&%ld&%ld",
							z->versionNumber,
							z->structSize,
							z->hashFunction,
							z->hash,
							z->salt,
							z->userID,
							z->userPW,
							z->computerID,
							z->exponent,
							z->days,
							z->requestType,
							z->programType,
							z->iteration,
							z->nextMsg);
			break;
		}
	case PRIMENET_RESULT_MESSAGE:		/* log a text message */
		{
			struct primenetResultMessage *z = (struct primenetResultMessage*) pkt;

			if ( validateID( z ) ) return PRIMENET_ERROR_INVALID_PARAM;

			armor( z->message );

			sprintf(args,"rm&%d&%d&%d&%u&%u&%s&%s&%s&%s",
							z->versionNumber,
							z->structSize,
							z->hashFunction,
							z->hash,
							z->salt,
							z->userID,
							z->userPW,
							z->computerID,
							z->message );
			break;
		}
	case PRIMENET_ASSIGNMENT_RESULT:	/* assignment result */
		
		{
			struct primenetAssignmentResult *z = (struct primenetAssignmentResult*) pkt;

			if ( validateID( z ) ) return PRIMENET_ERROR_INVALID_PARAM;

			if ( z->resultType == PRIMENET_RESULT_NOFACTOR )

				sprintf(args,"ar&%d&%d&%d&%u&%u&%s&%s&%s&%ld&%d&%04.1f",
							z->versionNumber,
							z->structSize,
							z->hashFunction,
							z->hash,
							z->salt,
							z->userID,
							z->userPW,
							z->computerID,
							z->exponent,
							z->resultType,
							z->resultInfo.how_far_factored);

			else if ( z->resultType == PRIMENET_RESULT_FACTOR ||
					  z->resultType == PRIMENET_RESULT_TEST ||
				      z->resultType == PRIMENET_RESULT_PRIME )

				sprintf(args,"ar&%d&%d&%d&%u&%u&%s&%s&%s&%ld&%d&%s",
							z->versionNumber,
							z->structSize,
							z->hashFunction,
							z->hash,
							z->salt,
							z->userID,
							z->userPW,
							z->computerID,
							z->exponent,
							z->resultType,
							z->resultInfo.factor);	/* all return szString[] from same offset */

			else if ( z->resultType == PRIMENET_RESULT_UNRESERVE )

				sprintf(args,"ar&%d&%d&%d&%u&%u&%s&%s&%s&%ld&%d",
							z->versionNumber,
							z->structSize,
							z->hashFunction,
							z->hash,
							z->salt,
							z->userID,
							z->userPW,
							z->computerID,
							z->exponent,
							z->resultType);

			else return PRIMENET_ERROR_INVALID_PARAM;

			break;
		}
	case PRIMENET_SET_COMPUTER_INFO:		/* set client's machine information */
		{
			struct primenetComputerInfo *z = (struct primenetComputerInfo*) pkt;

			if ( validateID( z ) ) return PRIMENET_ERROR_INVALID_PARAM;

			sprintf(args,"mi&%d&%d&%d&%u&%u&%s&%s&%s&%d&%d&%d",
							z->versionNumber,
							z->structSize,
							z->hashFunction,
							z->hash,
							z->salt,
							z->userID,
							z->userPW,
							z->computerID,
							z->cpu_type,
							z->speed,
							z->hours_per_day);
			break;
		}
	default:
		return PRIMENET_ERROR_INVALID_FIELD;	/* default error */
	}
	return 0;
}


/*////////////////////////////////////////////////////////////////////////////
//
// HTTP downloaded response page parsing procedures
//
/////////////////////////////////////////////////////////////////////////////*/

/* skip over the token name and point to the data string */

char* skip_token( char* s )
{
	while ( *s && *s != '=' ) s++;
	if ( *s == '=' ) s++;
	return s;
}


/* copy the data string up to the next '\r' delimiter character */

char* copy_value( char* buf, char* s )
{
	while ( *s && *s != '\r' ) *buf++ = *s++;
	if ( *s == '\r' ) s++;
	*buf = 0;
	return s;
}


/*
// parse_page: reads the server response page tokens and values
//             and converts these back into a PrimeNet 3.1 v4 packet
*/

int parse_page( char* buf, short operation, void* pkt )
{
	char *s, *stop;
	char lbuf[24];
	error_status_t res;

	/* get result code, which is always first */
	s = skip_token(buf);
	if ( !(*s) || *s < '0' || *s > '9' )
		return PRIMENET_ERROR_HTTP_BAD_PAGE;

	res = strtol( s, &s, 10 );

	/* quit now if garbage (for some reason) follows result code */
	if ( *s != '\r' ) return res;

	switch ( operation )
	{
			/* update completion, results message, assignment result */
			/* each have no data other than a result code */
	case PRIMENET_COMPLETION_DATE:
	case PRIMENET_ASSIGNMENT_RESULT:
	case PRIMENET_RESULT_MESSAGE:
	case PRIMENET_SET_COMPUTER_INFO:
			break;
	case PRIMENET_PING_SERVER_INFO:
		{
			/* ping returns server info */

			struct primenetPingServerInfo *z = (struct primenetPingServerInfo*) pkt;

			s = copy_value( lbuf, skip_token(s) );
			z->u.serverInfo.versionNumber = (short)strtol( lbuf, &stop, 10 );

			if ( z->u.serverInfo.versionNumber >= PRIMENET_PING_INFO_MESSAGE )
			{	/* v19+ message info ping packet */
				s = copy_value( lbuf, skip_token(s) );
				z->u.messageInfo.msgID = (short)strtol( lbuf, &stop, 10 );
				s = copy_value( lbuf, skip_token(s) );
				z->u.messageInfo.targetMap = strtoul( lbuf, &stop, 10 );
				s = copy_value( z->u.messageInfo.szMsgText, skip_token(s) );
#ifdef _DEBUG
				printf("versionNumber=%d\n", z->u.serverInfo.versionNumber);
				printf("msgID=%d\n", z->u.messageInfo.msgID);
				printf("targetMap=%08x\n", z->u.messageInfo.targetMap);
				printf("szMsgText=%s\n", z->u.messageInfo.szMsgText);
#endif
			}
			else
			{	/* normal server info ping packet */
				s = copy_value( z->u.serverInfo.buildID, skip_token(s) );
				s = copy_value( z->u.serverInfo.primenetServerName, skip_token(s) );
				s = copy_value( z->u.serverInfo.adminEmailAddr, skip_token(s) );
#ifdef _DEBUG
				printf("versionNumber=%d\n", z->u.serverInfo.versionNumber);
				printf("buildID=%s\n", z->u.serverInfo.buildID);
				printf("primenetServerName=%s\n", z->u.serverInfo.primenetServerName);
				printf("adminEmailAddr=%s\n", z->u.serverInfo.adminEmailAddr);
#endif
			}
			break;
		}
	case PRIMENET_MAINTAIN_USER_INFO:
		{
			/* update user info can return userID, userPW, computerID (usually what's sent) */

			struct primenetUserInfo *z = (struct primenetUserInfo*) pkt;

			s = copy_value( z->userID, skip_token(s) );
			s = copy_value( z->userPW, skip_token(s) );
			s = copy_value( z->computerID, skip_token(s) );

			if ( z->userID[0] == '.' ) z->userID[0] = 0;
			if ( z->userPW[0] == '.' ) z->userPW[0] = 0;
			if ( z->computerID[0] == '.' ) z->computerID[0] = 0;

#ifdef _DEBUG
			printf("userID=%s\n", z->userID);
			printf("userPW=%s\n", z->userPW);
			printf("computerID=%s\n", z->computerID);
#endif
			break;
		}
	case PRIMENET_GET_ASSIGNMENT:
		{
			/* get assignment returns requestType, exponent, how_far_factored */
	
			struct primenetGetAssignment *z = (struct primenetGetAssignment*) pkt;

			s = copy_value( lbuf, skip_token(s) );
			z->requestType = (short)strtol( lbuf, &stop, 10 );

			s = copy_value( lbuf, skip_token(s) );
			z->exponent = strtol( lbuf, &stop, 10 );

			s = copy_value( lbuf, skip_token(s) );
			z->how_far_factored = strtod( lbuf, &stop );

#ifdef _DEBUG
			printf("requestType=%d\n", z->requestType);
			printf("exponent=%ld\n", z->exponent);
			printf("how_factored=%4.1f\n", z->how_far_factored);
#endif
			break;
		}
	default:
		return PRIMENET_ERROR_HTTP_BAD_PAGE;
	}
	return res;
}


/*
// Primenet: main interface to Prime95.exe, accepts
//           and returns PrimeNet 3.0 packets
*/


int PrimeNet(short operation, void *pkt)
{
	error_status_t status = RPC_S_OK;
	char args[1024];		/* formatted arguments buffer */
	char pbuf[4096];		/* return page buffer */

	/* assemble GET/POST arguments */

	status = format_args( args, operation, pkt );
	if ( status != RPC_S_OK ) return status;

	/* send arguments, read back resulting page */

	status = pnHttpServer(pbuf, sizeof(pbuf), args);
	if ( status != RPC_S_OK ) return status;

	/* extract results from return page into packet */

	return parse_page( pbuf, operation, pkt );
}

/* EOF */

int (*PRIMENET)(short, void *) = PrimeNet;

/* Load the PrimeNet DLL, make sure an internet connection is active */

int LoadPrimeNet (void)
{
	/* Init stuff */
	/* Set PRIMENET procedure pointer */
	/* return false if not connected to internet */

	int lines = 0;
#ifndef AOUT
	FILE* fd;
	char buffer[4096];
#ifdef __EMX__
	char command[128];
	char szProxyHost[120], *con_host;
	char *colon;
	
	IniGetString(INIFILENAME, "ProxyHost", szProxyHost, 120, NULL);
	if (*szProxyHost) {
		if ((colon = strchr(szProxyHost, ':'))) {
			*colon = 0;
		}
		con_host = szProxyHost;
	} else {
		con_host = szSITE;
	}

	sprintf(command,"host %s",con_host);
#ifdef _DEBUG
	fprintf(stderr,"Command = %s\n",command);
#endif
	fd = popen(command,"r");
	if (fd != NULL) {
	  fgets(buffer, 199, fd);
#ifdef _DEBUG
	  fprintf(stderr,"Response = %s\n",buffer);
#endif
	  if (strncmp(buffer,"host:",5) != 0) {
	    fclose(fd);
	    return TRUE;
	  }
	  fclose(fd);
	}
#else
#ifdef __linux__
	/* Open file that will hopefully tell us if we are connected to */
	/* the Internet.  There are four possible settings for RouteRequired */
	/* 0:	Always return TRUE */
	/* 1:   Use old version 19 code */
	/* 2:   Use new code supplied by Matthew Ashton. */
	/* 99:	Default.  Use case 2 above but if cannot open /proc/net/route*/
	/*	then assume you are connected (we probably do not have read */
	/*	permission or this is a funny Linux setup). */
	{
	  int RtReq = IniGetInt (INIFILENAME, "RouteRequired", 99);
	  if (RtReq == 0) return (TRUE);
	  fd = fopen("/proc/net/route","r");
	  if (fd == NULL) return (RtReq == 99);
	/* We have a readable /proc/net/route file.  Use the new check */
	/* for an Internet connection written by Matthew Ashton. However, */
	/* we still support the old style check (just in case) by setting */
	/* RouteRequired to 1. */
	  if (RtReq >= 2) {
	    while (fgets(buffer, sizeof(buffer), fd)) {
	      int dest;
	      if(sscanf(buffer, "%*s %x", &dest) == 1 && dest == 0) {
		fclose (fd);
		return (TRUE);
	      }
	    }
	  }
	/* The old code for testing an Internet connection is below */
	  else {
	    fgets(buffer, 199, fd);
	    fgets(buffer, 199, fd);
	    while (!feof(fd)) {
	      if (strncmp(buffer, "lo", 2)) {
	        fclose(fd);
	        return TRUE;
	      }
	      fgets(buffer, 199, fd);
	    }
	  }
	  fclose(fd);
	}
#endif
#ifdef __FreeBSD__
	/* The /proc/net/route test is not really meaningful under FreeBSD */
	/* There doesn't seem to be any meaningful test to see whether the */
	/* computer is connected to the Internet at the time using a non- */
	/* invasive test (which wouldn't, say, activate diald or ppp or */
	/* something else */
	return TRUE;
#endif                /* __FreeBSD__ */
#endif
#endif
	OutputStr ("You are not connected to the Internet.\n");
	return FALSE;
}
