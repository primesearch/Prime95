/* Common strings */

#define MANUAL_QUIT	"You have elected to remove this computer from the Great Internet Mersenne Prime Search.  Other computers using this user ID will not be affected.  Please send the file results.txt to woltman@alum.mit.edu.\n\nAre you sure you want to do this?"
#define PRIMENET_QUIT	"You have elected to remove this computer from the Great Internet Mersenne Prime Search.  Other computers using this user ID will not be affected.\n\nPlease make sure your results have been successfully sent to the server (the program will be idle rather than looping trying to contact the server) before uninstalling the program. If in doubt, you can send the results.txt file to woltman@alum.mit.edu.\n\nYou can either complete your current assignment or you can quit GIMPS immediately.  Do you wish to complete your current work assignments before quitting?"
#define PING_ERROR	"Unable to get version information from PrimeNet server."
#define MSG_MEMORY  "You have left the available memory fields at 8 megabytes.  You can increase your chances of finding a Mersenne prime slightly if you let prime95 occasionally use more memory.  The readme.txt file has more information.  Do you want to let prime95 use more memory?"
#define MSG_SPEED  "You have entered a CPU speed that may not be correct.  An incorrect CPU speed entry will result in inaccurate timings.  Are you sure this is the correct CPU speed value?"

/* Global variables */

extern int STARTUP_IN_PROGRESS;	/* TRUE if we should call primeContinue */
				/* after asking for user information */

/* Common routines */

void sanitizeString (char *);
void rangeStatusMessage (char *);
