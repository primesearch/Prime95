/* Common strings */

#define MANUAL_QUIT	"You have elected to remove this computer from the Great Internet Mersenne Prime Search.  Other computers using this user ID will not be affected.  Please send the file results.txt to woltman@alum.mit.edu.\n\nAre you sure you want to do this?"
#define PRIMENET_QUIT	"You have elected to remove this computer from the Great Internet Mersenne Prime Search.  Other computers using this user ID will not be affected.  All unfinished work will be returned to the PrimeNet server.  Make sure you are connected to the Internet before choosing \"Yes\".\n\nAre you sure you want to do this?"
#define PING_ERROR	"Unable to get version information from PrimeNet server."
#define MSG_MEMORY  "You have left the available memory fields at 8 megabytes.  You can increase your chances of finding a Mersenne prime slightly if you let prime95 occasionally use more memory.  The readme.txt file has more information.  Do you want to let prime95 use more memory?"

/* Global variables */

extern int STARTUP_IN_PROGRESS;	/* TRUE if we should call primeContinue */
				/* after asking for user information */

/* Common routines */

void sanitizeString (char *);
void rangeStatusMessage (char *);
