/* Include file for PRIMENET.DLL interface */

/* user authentication string limits */

#define PRIMENET_USER_ID_LENGTH	14
#define PRIMENET_USER_PW_LENGTH	8
#define PRIMENET_COMPUTER_ID_LENGTH 12

/* Each program can use a different private hash function */

#define PRIMENET_HASHFUNC_WOLTMAN_1	4

/* The current version number of the the PRIMENET.DLL API */

#define	PRIMENET_VERSION	4

/* There is only one routine in PRIMENET.DLL.  This routine, */
/* "PrimeNet", takes two arguments one is an operation_type */
/* defined below, and the other argument is a pointer to a */
/* structure.  A different structure is declared for each */
/* operation type.  All strings in structures are zero-terminated. */

/* The PrimeNet routine is declared as such: */
/*	int __stdcall *PrimeNet (short operation, void *pkt); */
/* Prime95 loads the routine as follows: */
/*	HPROC_PRIMENET =					*/
/*		(int (__stdcall *)(short, void *))		*/
/*			GetProcAddress (HLIB, "PrimeNet");	*/

#define PRIMENET_MAINTAIN_USER_INFO	0
#define PRIMENET_GET_ASSIGNMENT		1
#define PRIMENET_COMPLETION_DATE	2
#define PRIMENET_RESULT_MESSAGE		3
#define PRIMENET_ASSIGNMENT_RESULT	4
#define PRIMENET_PING_SERVER_INFO	5
#define PRIMENET_SET_COMPUTER_INFO	6


/* This structure is passed for the PING_SERVER_INFO call */
/* The server replies with server information or a message */
/* that the server is broadcasting to all clients.  The server can */
/* target specific client versions and OSes.  The targetMap field */
/* combines the VERSION_BIT from commonc.h to target specific client */
/* versions and the PORT bit to target spefic OSes. */

#define PRIMENET_PING_INFO_MESSAGE	0x7F00

struct primenetPingServerInfo {
	union {				/* size is 202 */
		struct serverInfo {	/* if versionNumber < 0x7F00 */
			short	versionNumber;
			char	buildID[40];
			char	primenetServerName[80];
			char	adminEmailAddr[80];
		} serverInfo;
		struct messageInfo {	/* if versionNumber >= 0x7F00 */
			short	versionNumber;
			short	msgID;
			long	targetMap;
			char	szMsgText[194];
		} messageInfo;
	} u;
};

#define PRIMENET_USER_OPTION_SENDEMAIL	0x01
#define PRIMENET_USER_OPTION_TEAMACCT	0x02

/* This structure is passed for the MAINTAIN_USER_INFO call */

struct primenetUserInfo {
	short	versionNumber;
	short	structSize;		/* Size of this structure */
	short	hashFunction;		/* Each program type uses a */
					/* different hash function */
	unsigned short hash;		/* A hash value of all following */
					/* fields (for security) */
	unsigned short salt;		/* A random value to make cracking */
					/* security codes even harder */
	char	userID[PRIMENET_USER_ID_LENGTH+1];
	char	userPW[PRIMENET_USER_PW_LENGTH+1];
	char	computerID[PRIMENET_COMPUTER_ID_LENGTH+1];
	char	userName[80];
	char	userEmailAddr[80];
	char	oldUserID[PRIMENET_USER_ID_LENGTH+1];
					/* Only used when changing from */
	char	oldUserPW[PRIMENET_USER_PW_LENGTH+1];
					/* one userid to another userid */
	char	bUserOptions;		/* TRUE if server should e-mail */
					/* exponent about to expire notices */
};


/* This structure is passed for the GET_ASSIGNMENT call */

#define PRIMENET_ASSIGN_FACTOR			0x0001
#define PRIMENET_ASSIGN_TEST			0x0002
#define PRIMENET_ASSIGN_DBLCHK			0x0004
#define PRIMENET_ASSIGN_DBLCHK_SUSPICIOUS	0x0008
#define PRIMENET_ASSIGN_BIGONES			0x0010
#define PRIMENET_ASSIGN_PFACTOR			0x0020

#define PRIMENET_PROGRAM_OLDWOLTMAN		0x0001
#define PRIMENET_PROGRAM_WOLTMAN		0x0002
#define PRIMENET_PROGRAM_MERS			0x0004
#define PRIMENET_PROGRAM_MACLUCAS		0x0008
#define PRIMENET_PROGRAM_OTHER			0x0010

struct primenetGetAssignment {
	short	versionNumber;
	short	structSize;		/* Size of this structure */
	short	hashFunction;		/* Each program type uses a */
					/* different hash function */
	unsigned short hash;		/* A hash value of all following */
					/* fields (for security) */
	unsigned short salt;		/* A random value to make cracking */
					/* security codes even harder */
	char	userID[PRIMENET_USER_ID_LENGTH+1];
	char	userPW[PRIMENET_USER_PW_LENGTH+1];
	char	computerID[PRIMENET_COMPUTER_ID_LENGTH+1];
	short	requestType;		/* One or more of the request */
					/* types above. */
	short	programType;		/* One of the program types.  Only */
					/* used by server when requesting */
					/* a double-check assignment */
	unsigned long exponent;		/* Returned by the server */
	double	how_far_factored;	/* Log base 2 of highest trial */
					/* factor tested. */
};


/* This structure is passed for the COMPLETION_DATE call */
/* This must be called after every getAssignment call.  It may */
/* also be called whenever prime95 feels like it. */

struct primenetCompletionDate {
	short	versionNumber;
	short	structSize;		/* Size of this structure */
	short	hashFunction;		/* Each program type uses a */
					/* different hash function */
	unsigned short hash;		/* A hash value of all following */
					/* fields (for security) */
	unsigned short salt;		/* A random value to make cracking */
					/* security codes even harder */
	char	userID[PRIMENET_USER_ID_LENGTH+1];
	char	userPW[PRIMENET_USER_PW_LENGTH+1];
	char	computerID[PRIMENET_COMPUTER_ID_LENGTH+1];
	unsigned long exponent;		/* The checked out exponent */
	unsigned long days;		/* Number of days to completion */
	short	requestType;		/* Server wants these for */
	short	programType;		/* unexpected checkins */
	unsigned long iteration;	/* Current LL iter or factoring pass */
	unsigned long nextMsg;		/* Days until next expected */
					/* completion date msg will be sent */
};


/* This structure is passed for the RESULT_MESSAGE call */

struct primenetResultMessage {
	short	versionNumber;
	short	structSize;		/* Size of this structure */
	short	hashFunction;		/* Each program type uses a */
					/* different hash function */
	unsigned short hash;		/* A hash value of all following */
					/* fields (for security) */
	unsigned short salt;		/* A random value to make cracking */
					/* security codes even harder */
	char	userID[PRIMENET_USER_ID_LENGTH+1];
	char	userPW[PRIMENET_USER_PW_LENGTH+1];
	char	computerID[PRIMENET_COMPUTER_ID_LENGTH+1];
	char	message[200];
	char	pad1[1];		/* Make struct size even */
};


/* This structure is passed for the ASSIGNMENT_RESULT call */

/* When the server receives a "no factor" message, it updates */
/* its database and frees the reservation if the exponent was reserved */
/* for factoring only.  When the server receives a Lucas-Lehmer test */
/* completed message, it removes the exponent from the database */
/* The exponent may later be added back as needing a double-check. */
/* When the server receives a factor message, it removes the exponent */
/* from the database.  When the server receives a unreserve message, */
/* it frees the exponent for others to reserve and test. */

#define PRIMENET_RESULT_NOFACTOR	0
#define PRIMENET_RESULT_TEST		1
#define PRIMENET_RESULT_FACTOR		2
#define PRIMENET_RESULT_PRIME		3
#define PRIMENET_RESULT_UNRESERVE	4

struct primenetAssignmentResult {
	short	versionNumber;
	short	structSize;		/* Size of this structure */
	short	hashFunction;		/* Each program type uses a */
					/* different hash function */
	unsigned short hash;		/* A hash value of all following */
					/* fields (for security) */
	unsigned short salt;		/* A random value to make cracking */
					/* security codes even harder */
	char	userID[PRIMENET_USER_ID_LENGTH+1];
	char	userPW[PRIMENET_USER_PW_LENGTH+1];
	char	computerID[PRIMENET_COMPUTER_ID_LENGTH+1];
	char	pad1[1];		/* For proper alignment */
	unsigned long exponent;		/* Exponent that was tested */
	short	resultType;		/* One of the above result types */
	char	pad2[2];		/* For proper alignment */
	union primenetResultInfo {
		double	how_far_factored;
		char	residue[17];
		char	factor[32];
	} resultInfo;
};

/* Cpu types */

#define PRIMENET_CPU_CYRIX		3
#define PRIMENET_CPU_486		4
#define PRIMENET_CPU_PENTIUM		5
#define PRIMENET_CPU_PPRO_OR_PII	6
#define PRIMENET_CPU_K6			7

/* This structure is passed for the SET_COMPUTER_INFO call */

struct primenetComputerInfo {
	short	versionNumber;
	short	structSize;		/* Size of this structure */
	short	hashFunction;		/* Each program type uses a */
					/* different hash function */
	unsigned short hash;		/* A hash value of all following */
					/* fields (for security) */
	unsigned short salt;		/* A random value to make cracking */
					/* security codes even harder */
	char	userID[PRIMENET_USER_ID_LENGTH+1];
	char	userPW[PRIMENET_USER_PW_LENGTH+1];
	char	computerID[PRIMENET_COMPUTER_ID_LENGTH+1];
	short	cpu_type;		/* enumerated CPU model */
	short	speed;			/* speed in Mhz */
	short	hours_per_day;		/* 0 to 24 */
};



/* Error codes returned to Prime95 */

#define PRIMENET_ERROR_OK		0	/* no error */
#define PRIMENET_ERROR_SERVER_UNSPEC	1	/* unspecified server error */
#define PRIMENET_ERROR_UNASSIGNED	3	/* exponent is not assigned */
#define PRIMENET_ERROR_ACCESS_DENIED	5	/* unauthorized action */
#define PRIMENET_ERROR_NO_ASSIGNMENT	7	/* no assignment matches */
#define PRIMENET_ERROR_ALREADY_TESTED	11	/* result already in server */
#define PRIMENET_ERROR_INVALID_FIELD	13	/* DLL got bad data in packet*/
#define PRIMENET_ERROR_NOT_PERMITTED	17	/* exp owned by someone else */
#define PRIMENET_ERROR_SERVER_BUSY	23	/* server busy */
#define PRIMENET_ERROR_HTTP_BAD_PAGE	29	/* HTTP response not PrimeNet*/
						/* formatted */
#define PRIMENET_ERROR_HTTP_BAD_PACKET	31	/* HTTP packet received is */
						/* invalid */
#define PRIMENET_ERROR_INVALID_PARAM	87	/* invalid parameter */
#define PRIMENET_ERROR_CONNECT_FAILED	2250	/* server unavailable */

#define PRIMENET_ERROR_MODEM_OFF	12345	/* modem offline */
#define PRIMENET_ERROR_BLACKOUT		12346	/* blackout! */

#define PRIMENET_ERROR_KILL_CLIENT	31000	/* client not supported by */
						/* primenet anymore */
#define PRIMENET_ERROR_SLEEP_CLIENT_1	32001	/* client should not contact */
#define PRIMENET_ERROR_SLEEP_CLIENT_2	32002	/* primenet server for next */
#define PRIMENET_ERROR_SLEEP_CLIENT_200	32200	/* 1,2,...,200 hours. */

/* Old packets */

struct primenet2ResultMessage {
	short	versionNumber;
	char	userID[PRIMENET_USER_ID_LENGTH+1];
	char	userPW[PRIMENET_USER_PW_LENGTH+1];
	char	message[200];
};

struct primenet2AssignmentResult {
	short	versionNumber;
	char	userID[PRIMENET_USER_ID_LENGTH+1];
	char	userPW[PRIMENET_USER_PW_LENGTH+1];
	unsigned long exponent;
	short	resultType;
	union primenet2ResultInfo {
		double	how_far_factored;
		char	residue[17];
		char	factor[32];
	} resultInfo;
};
