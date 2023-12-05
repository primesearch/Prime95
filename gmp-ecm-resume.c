





/* Functions for reading a writing resume file lines.

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2010, 2011, 2012
Paul Zimmermann, Alexander Kruppa and Cyril Bouvier.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#if !defined (_MSC_VER)
#include <unistd.h>
#endif
#include <gmp.h>
//#include "ecm.h"
//#include "ecm-ecm.h"

#define EC_W_NBUFS 9 /* for Hessian form */

#define ECM_EC_TYPE_MONTGOMERY           1
#define ECM_EC_TYPE_WEIERSTRASS          2
#define ECM_EC_TYPE_HESSIAN              3
#define ECM_EC_TYPE_WEIERSTRASS_COMPLETE 4

/* Different parametrizations used in stage 1 of ECM */
#define ECM_PARAM_DEFAULT -1
#define ECM_PARAM_SUYAMA 0
#define ECM_PARAM_BATCH_SQUARE 1
#define ECM_PARAM_BATCH_2 2
#define ECM_PARAM_BATCH_32BITS_D 3
/* we keep 4 as spare */
#define ECM_PARAM_WEIERSTRASS 5
#define ECM_PARAM_HESSIAN 6
#define ECM_PARAM_TORSION 7

/* different methods implemented */
#define ECM_ECM 0
#define ECM_PM1 1
#define ECM_PP1 2

/* The checksum for savefile is the product of all mandatory fields, modulo
   the greatest prime below 2^32 */
#define CHKSUMMOD 4294967291U

/* Apple uses '\r' for newlines */
#define IS_NEWLINE(c) (((c) == '\n') || ((c) == '\r'))

/* Structure for candidate usage.  This is much more powerful than using a
   simple mpz_t to hold the candidate.  This structure also houses the 
   expression (in raw form), and will modify the expression as factors 
   are found (if in looping modes).  Also, since we are warehousing all
   of the data associated with the candidate, we also store whether the
   candidate is PRP here (so testing will cease), along with the length
   of the candidate.  As each factor is found, the candidate will also
   have the factor removed from it */
typedef struct
{
#if defined (CANDI_DEBUG)
  unsigned long magic;	/* used for debugging purposes while writing this code */
#endif
  char *cpExpr;		/* if non-NULL, then this is a "simpler" expression than the 
			   decimal output of n */
  mpz_t n;		/* the cofactor candidate currently being used to find factors from */
  unsigned ndigits;	/* the number of digits (decimal) in n */
  unsigned nexprlen;	/* strlen of expression, 0 if there is NO expression */
  int isPrp;		/* usually 0, but turns 1 if factor found, and the cofactor is PRP, 
			   OR if the original candidate was PRP and the user asked to prp check */
} mpcandi_t;

typedef struct
{
  int type;              
  int law;
  mpz_t a4;               /* for MONTGOMERY: b*y^2=x^3+A*x^2+x 
			      for WEIERSTRASS: y^2=x^3+A*x+B
			      for HESSIAN: U^3+V^3+W^3=3*A*U*V*W */
  mpz_t a1, a3, a2, a6;  /* for complete WEIERSTRASS */
  mpz_t buf[EC_W_NBUFS]; /* used in the addition laws */
  int disc;                /* in case E is known to have CM by Q(sqrt(disc)) */
  mpz_t sq[10];          /* for CM curves, we might have squareroots */
} __ell_curve_struct;
typedef __ell_curve_struct ell_curve_t[1];

typedef struct
{
  int method;     /* factorization method, default is ecm */
  mpz_t x, y;        /* starting point (if non zero) */
  int param;      /* (ECM only) What parametrization do we used */ 
  mpz_t sigma;    /* (ECM only) The parameter for the parametrization */ 
                      /* May contains A */
  int sigma_is_A; /* if  1, 'parameter' contains A (Montgomery form),
		     if  0, 'parameter' contains sigma (Montgomery form),
		     if -1, 'parameter' contains A, and the input curve is in
		     Weierstrass form y^2 = x^3 + A*x + B, with y in 'go'. */
  __ell_curve_struct *E;   /* the curve, particularly useful for CM ones */
  mpz_t go;       /* initial group order to preload (if NULL: do nothing),
		     or y for Weierstrass form if sigma_is_A = -1. */
  double B1done;  /* step 1 was already done up to B1done */
  mpz_t B2min;    /* lower bound for stage 2 (default is B1) */
  mpz_t B2;       /* step 2 bound (chosen automatically if < 0.0) */
  unsigned long k;/* number of blocks in stage 2 */
  int S;          /* degree of the Brent-Suyama's extension for stage 2 */
  int repr;       /* representation for modular arithmetic: ECM_MOD_MPZ=mpz,         
		     ECM_MOD_MODMULN=modmuln (Montgomery's quadratic multiplication),
		     ECM_MOD_REDC=redc (Montgomery's subquadratic multiplication),
		     ECM_MOD_GWNUM=Woltman's gwnum routines (tbd),
		     > 16 : special base-2 representation        
		     MOD_DEFAULT: automatic choice */
  int nobase2step2; /* disable special base-2 code in ecm stage 2 only */
  int verbose;    /* verbosity level: 0 no output, 1 normal output,   
		     2 diagnostic output */
  FILE *os;       /* output stream (for verbose messages) */
  FILE *es;       /* error  stream (for error   messages) */
  char *chkfilename; /* Filename to write stage 1 checkpoints to */
  char *TreeFilename; /* Base filename for storing product tree of F */
  double maxmem;  /* Maximal amount of memory to use in stage 2, in bytes.
                     0. means no limit (optimise only for speed) */
  double stage1time; /* Time to add for estimating expected time to find fac.*/
  gmp_randstate_t rng; /* State of random number generator */
  int use_ntt;     /* set to 1 to use ntt poly code in stage 2 */
  int (*stop_asap) (void); /* Pointer to function, if it returns 0, contine 
                      normally, otherwise exit asap. May be NULL */
  /* The batch mode is used for stage 1 when param=1 or param=2)*/
  mpz_t batch_s;   /* s is the product of primes up to B1 for batch mode */
  double batch_last_B1_used; /* Last B1 used in batch mode. Used to avoid */
                             /*  computing s when B1 = batch_last_B1_used */
  int gpu;  /* do we use the GPU for stage 1. */
            /* If different from 0, the GPU is used */
            /* Else, the parameters beginning by gpu_* have no meaning */
  int gpu_device; /* Which device do we use */
  int gpu_device_init; /* Is the device initialized?*/
  unsigned int gpu_number_of_curves; 
  double gw_k;         /* use for gwnum stage 1 if input has form k*b^n+c */
  unsigned long gw_b;  /* use for gwnum stage 1 if input has form k*b^n+c */
  unsigned long gw_n;  /* use for gwnum stage 1 if input has form k*b^n+c */
  signed long gw_c;    /* use for gwnum stage 1 if input has form k*b^n+c */
} __ecm_param_struct;
typedef __ecm_param_struct ecm_params[1];
typedef __ecm_param_struct *ecm_params_ptr;

//#ifdef HAVE_FCNTL_H
#include <fcntl.h>
//#endif

#if defined (_MSC_VER) || defined (__MINGW32__)
/* needed to declare GetComputerName() for write_resumefile_line() */
#include <windows.h>
#endif



/* Tries to read a number from a line from fd and stores it in r.
   Keeps reading lines until a number is found. Lines beginning with "#"
     are skipped.
   Returns 1 if a number was successfully read, 0 if no number can be read
     (i.e. at EOF)
   Function is now simpler.  Much of the logic (other than skipping # lines
     is now contained within eval() function.
*/

int
read_number (mpcandi_t *n, FILE *fd, int primetest)
{
  int c;

new_line:
  c = fgetc (fd);

  /* Skip comment lines beginning with '#' */
  if (c == '#')
    {
      do
        c = fgetc (fd);
      while (c != EOF && !IS_NEWLINE(c));
      if (IS_NEWLINE(c))
        goto new_line;
    }

  if (c == EOF)
    return 0;

  ungetc (c, fd);
//GW:
//  if (!eval (n, fd, primetest))
//    goto new_line;

#if 0
  /*  Code to test out eval_str function, which "appears" to work correctly. */
  {
    /* warning!! Line is pretty small, but since this is just testing code, we
       can easily control the input for this test.  This code should NEVER be
       compiled into released build, its only for testing of eval_str() */
    char Line[500], *cp;
    fgets (Line, sizeof(Line), fd);

    if (!eval_str (n, Line, primetest, &cp))
      goto new_line;
    fprintf (stderr, "\nLine is at %X cp is at %X\n", Line, cp);
  }
#endif

#if defined (DEBUG_EVALUATOR)
  if (n->cpExpr)
    fprintf (stderr, "%s\n", n->cpExpr);
  mpz_out_str (stderr, 10, n->n);
  fprintf (stderr, "\n");
#endif

  return 1;
}




/* Reads a string of characters from fd while they match the string s.
   Returns the number of matching characters that were read. 
*/

static int 
facceptstr (FILE *fd, char *s)
{
  int c;
  unsigned i = 0;
  
  while (s[i] != 0 && (c = fgetc (fd)) != EOF)
    {
      if (c != s[i++])
        {
          ungetc (c, fd);
          return i-1;
        }
    }
  
  return i;
}

/* Accepts "\n" or "\r\n" or "\r". 
   Returns 1 if any of the three was read, 0 otherwise */

static int 
facceptnl (FILE *fd)
{
  int c, r = 0;

  c = fgetc (fd);
  if (c == '\r')
    {
      c = fgetc (fd);
      r = 1;
    }

  if (c == '\n')
    r = 1;
  else if (c != EOF)
    ungetc (c, fd);

  return r;
}

/* Reads a string from fd until the character "delim" or newline is seen, or 
   "len" characters have been written to "s" (including terminating null), 
   or EOF is reached. The "delim" and newline characters are left on the 
   stream.
   If s is NULL, characters are read from fd but not written anywhere.
   Returns the number of characters read.
*/

static int 
freadstrn (FILE *fd, char *s, char delim, unsigned int len)
{
  unsigned int i = 0;
  int c;
  
  while (i + 1 < len && (c = fgetc (fd)) != EOF)
    if (c == delim || IS_NEWLINE(c))
      {
        ungetc (c, fd);
        break;
      }
    else
      if (s != NULL)
        s[i++] = (char) c;
  
  if (i < len && s != NULL)
    s[i++] = 0;
  
  return i;
}

/* Reads an assignment from a save file. Return 1 if an assignment was
   successfully read, 0 if there are no more lines to read (at EOF) 
*/

int 
read_resumefile_line (int *method, mpz_t x, mpz_t y, mpcandi_t *n, 
		      mpz_t sigma, mpz_t A,
		      mpz_t x0, mpz_t y0, int *Etype, int *param, 
		      double *b1, char *program, char *who, char *rtime, 
		      char *comment, FILE *fd)
{
  int a, have_method, have_x, have_y, have_z, have_n, have_sigma, have_a, 
      have_b1, have_checksum, have_qx;
  unsigned int saved_checksum;
  char tag[16];
  mpz_t z;
  
  while (!feof (fd))
    {
      /* Ignore empty lines */
      if (facceptnl (fd))
        continue;
      
      /* Ignore lines beginning with '#'*/
      if (facceptstr (fd, "#"))
        {
          while (!facceptnl (fd) && !feof (fd))
            fgetc (fd);
          continue;
        }
      
      if (feof (fd))
        break;
      
      have_method = have_x = have_y = have_z = have_n = have_sigma = have_a = 
                    have_b1 = have_qx = have_checksum = 0;

      /* For compatibility reason, param = ECM_PARAM_SUYAMA by default */
      *param = ECM_PARAM_SUYAMA;

      /* default and compatibility reasons */
      *Etype = ECM_EC_TYPE_MONTGOMERY;

      /* Set optional fields to zero */
      mpz_set_ui (sigma, 0);
      mpz_set_ui (A, 0);
      if (program != NULL)
        program[0] = 0;
      if (who != NULL)
        who[0] = 0;
      if (rtime != NULL)
        rtime[0] = 0;
      if (comment != NULL)
        comment[0] = 0;

      while (!facceptnl (fd) && !feof (fd))
        {
          freadstrn (fd, tag, '=', 16);
          
          if (!facceptstr (fd, "="))
            {
              printf ("Resume warning, skipping line with no '=' after: %s\n", tag);
              goto error;
            }
          
          if (strcmp (tag, "METHOD") == 0)
            {
              if (facceptstr (fd, "ECM") == 3)
                *method = ECM_ECM;
              else if (facceptstr (fd, "P"))
                {
                  a = facceptstr (fd, "-1");
                  if (a == 2)
                    *method = ECM_PM1;
                  else if (a == 0 && facceptstr (fd, "+1") == 2)
                    *method = ECM_PP1;
                  else
                    goto error;
                }
              else
                goto error;

              have_method = 1;
            }
          else if (strcmp (tag, "X") == 0)
            {
              mpz_inp_str (x, fd, 0);
              have_x = 1;
            }
          else if (strcmp (tag, "Y") == 0)
            {
              mpz_inp_str (y, fd, 0);
              have_y = 1;
            }
          else if (strcmp (tag, "Z") == 0)
            {
              mpz_init (z);
              mpz_inp_str (z, fd, 0);
              have_z = 1;
            }
          else if (strcmp (tag, "QX") == 0)
            {
              mpz_inp_str (x, fd, 0);
              have_qx = 1;
            }
          else if (strcmp (tag, "X0") == 0)
            {
              mpz_inp_str (x0, fd, 0);
            }
          else if (strcmp (tag, "Y0") == 0)
            {
              mpz_inp_str (y0, fd, 0);
            }
          else if (strcmp (tag, "CHECKSUM") == 0)
            {
              if (fscanf (fd, "%u", &saved_checksum) != 1)
                goto error;
              have_checksum = 1;
            }
          else if (strcmp (tag, "COMMENT") == 0)
            {
              freadstrn (fd, comment, ';', 255);
            }
          else if (strcmp (tag, "N") == 0)
            {
              /*mpz_inp_str (n, fd, 0);*/
	      /* we want to "maintain" any expressions, which were possibly stored in the file for N */
              have_n = read_number (n, fd, 0);
            }
          else if (strcmp (tag, "SIGMA") == 0)
            {
              mpz_inp_str (sigma, fd, 0);
              have_sigma = 1;
            }
          else if (strcmp (tag, "PARAM") == 0)
            {
              if (fscanf (fd, "%d", param) != 1)
                goto error;
            }
          else if (strcmp (tag, "ETYPE") == 0)
            {
              if (fscanf (fd, "%d", Etype) != 1)
                goto error;
            }
          else if (strcmp (tag, "A") == 0)
            {
              mpz_inp_str (A, fd, 0);
              have_a = 1;
            }
          else if (strcmp (tag, "B1") == 0)
            {
              if (fscanf (fd, "%lf", b1) != 1)
                goto error;
              have_b1 = 1;
            }
          else if (strcmp (tag, "PROGRAM") == 0)
            {
              freadstrn (fd, program, ';', 255);
            }
          else if (strcmp (tag, "WHO") == 0)
            {
              freadstrn (fd, who, ';', 255);
            }
          else if (strcmp (tag, "TIME") == 0)
            {
              freadstrn (fd, rtime, ';', 255);
            }
          else /* Not a tag we know about */
            {
              printf ("Save file line has unknown tag: %s\n", tag);
              goto error;
            }
         
          /* Prime95 lines have no semicolon after SIGMA */
          if (!facceptstr (fd, ";") && ! (have_qx && have_n && have_sigma))
            {
              printf ("%s field not followed by semicolon\n", tag);
              goto error;
            }
          
          while (facceptstr (fd, " "));
        }
      
      /* Finished reading tags */
      
      /* Handle Prime95 v22 lines. These have no METHOD=ECM field and
         QX= instead of X= */
      
      if (have_qx)
        {
          if (have_n && have_sigma)
            {
              *method = ECM_ECM;
              /* *b1 = 1.0; */
              strcpy (program, "Prime95");
              mpz_mod (x, x, n->n);
              return 1;
            }
          goto error;
        }

#ifdef DEBUG
      if (*method != ECM_ECM && (have_sigma || have_a || have_z))
        {
          int count = have_sigma + have_a + have_z;
          printf ("Warning: Save file line has");
          if (have_sigma)
            {
              printf (" SIGMA");
              mpz_set_ui (sigma, 0);
              if (--count > 1)
                printf (",");
              else if (count > 0)
                printf (" and");
            }
          if (have_a)
            {
              printf (" A");
              mpz_set_ui (A, 0);
              if (--count > 0)
                printf (" and");
            }
          if (have_z)
            {
              printf (" Z");
              mpz_clear (z);
              have_z = 0;
            }
          printf (" value for method other than ECM.\n");
        }
#endif
      
      if (!have_method || !have_x || !have_n || !have_b1 ||
          (*method == ECM_ECM && !have_sigma && !have_a))
        {
          fprintf (stderr, "Save file line lacks fields\n");
          continue;
        }

      if (have_checksum)
        {
          mpz_t checksum;
          
          mpz_init (checksum);
          mpz_set_d (checksum, *b1);
          if (have_sigma)
            mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (sigma, CHKSUMMOD));
          if (have_a)
            mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (A, CHKSUMMOD));
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (n->n, CHKSUMMOD));
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (x, CHKSUMMOD));
          if (have_z)
            mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (z, CHKSUMMOD));
          mpz_mul_ui (checksum, checksum, (*param+1)%CHKSUMMOD);
          if (mpz_fdiv_ui (checksum, CHKSUMMOD) != saved_checksum)
            {
              fprintf (stderr, "Resume file line has bad checksum %u, expected %lu\n", 
                       saved_checksum, mpz_fdiv_ui (checksum, CHKSUMMOD));
              mpz_clear (checksum);
              continue;
            }
          mpz_clear (checksum);
        }

      mpz_mod (x, x, n->n);
      if (have_y)
        mpz_mod(y, y, n->n);
      if (have_z)	/* Must normalize */
        {
          if (!mpz_invert (z, z, n->n)) /* Factor found? */
            {
              /* Oh great. What do we do with it now? */
              /* mpres_gcd (f, z, n); */
              printf ("Oops, factor found while reading from save file.\n");
            }
          mpz_mul (z, z, x);
          mpz_mod (x, z, n->n);
          mpz_clear (z);
        }

      return 1;
      
error:
      /* This can occur when reading Prime95 resume files,
         or files that have comment lines in them,
         or files that have a problem with the save line */
      /* In case of error, read rest of line and try next line */
      while (!facceptnl (fd) && !feof (fd))
        fgetc (fd);
    }
    
    /* We hit EOF without reading a proper save line */
    return 0;
}


/* Append a residue in file. */
static void  
write_resumefile_line (FILE *file, int method, double B1, mpz_t sigma, 
                       int sigma_is_A, int Etype, int param, mpz_t x, mpz_t y,
		       mpcandi_t *n, mpz_t x0, mpz_t y0, const char *comment)
{
  mpz_t checksum;
  time_t t;
  char text[256];
  char *uname, mname[32];

  mpz_init (checksum);
  mpz_set_d (checksum, B1);
  fprintf (file, "METHOD=");
  if (method == ECM_PM1)
    fprintf (file, "P-1");
  else if (method == ECM_PP1)
    fprintf (file, "P+1");
  else 
    {
      fprintf (file, "ECM");
      if (sigma_is_A == 0)
        {
          if (param != ECM_PARAM_DEFAULT)
            fprintf (file, "; PARAM=%d", param);

          fprintf (file, "; SIGMA=");
        }
      else
        fprintf (file, "; ETYPE=%d; A=", Etype);
          
        mpz_out_str (file, 10, sigma);
        mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (sigma, CHKSUMMOD));
        if (param != ECM_PARAM_DEFAULT)
            mpz_mul_ui (checksum, checksum, (param+1)%CHKSUMMOD);
    }
  
  fprintf (file, "; B1=%.0f; N=", B1);
  if (n->cpExpr)
    fprintf(file, "%s", n->cpExpr);
  else
    mpz_out_str (file, 10, n->n);
  fprintf (file, "; X=0x");
  mpz_out_str (file, 16, x);
  mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (n->n, CHKSUMMOD));
  mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (x, CHKSUMMOD));
  fprintf (file, "; CHECKSUM=%lu; PROGRAM=GMP-ECM %s;",
           mpz_fdiv_ui (checksum, CHKSUMMOD), VERSION);
  mpz_clear (checksum);
  if (y != NULL)
    {
      fprintf (file, " Y=0x");
      mpz_out_str (file, 16, y);
      fprintf (file, ";");
    }
  
  if (x0 != NULL)
    {
      fprintf (file, " X0=0x");
      mpz_out_str (file, 16, x0);
      fprintf (file, ";");
    }
  if (y0 != NULL)
    {
      fprintf (file, " Y0=0x");
      mpz_out_str (file, 16, y0);
      fprintf (file, ";");
    }
  
  /* Try to get the users and his machines name */
  /* TODO: how to make portable? */
  uname = getenv ("LOGNAME");
  if (uname == NULL)
    uname = getenv ("USERNAME");
  if (uname == NULL)
    uname = "";
  
#if defined (_MSC_VER) || defined (__MINGW32__)
  /* dummy block, so that the vars needed here don't need to
    "spill" over to the rest of the function. */
  {
    DWORD size;
    /* x86_64-w64-mingw32-gcc (GCC) 4.8.0 20121031 (experimental) has infinite
       for loop below with -O2, volatile seems to fix it */
    volatile size_t i;
    TCHAR T[MAX_COMPUTERNAME_LENGTH+2];
    size=MAX_COMPUTERNAME_LENGTH+1;
    if (!GetComputerName(T, &size))
      strcpy(mname, "localPC");
    else
      {
        if (size > sizeof(mname) - 1)
          size = sizeof(mname) - 1;

        for (i = 0; i < size; ++i)
          mname[i] = T[i];
        mname[i] = 0;
      }
  }
#else
  if (gethostname (mname, 32) != 0)
    mname[0] = 0;
  mname[31] = 0; /* gethostname() may omit trailing 0 if hostname >31 chars */
#endif
  
  if (uname[0] != 0 || mname[0] != 0)
    {
      fprintf (file, " WHO=%.233s@%.32s;", uname, mname);
    }

  if (comment[0] != 0)
    fprintf (file, " COMMENT=%.255s;", comment);
  
  t = time (NULL);
  strncpy (text, ctime (&t), 255);
  text[255] = 0;
  text[strlen (text) - 1] = 0; /* Remove newline */
  fprintf (file, " TIME=%s;", text);
  fprintf (file, "\n");
  fflush (file);
}

/* Call write_resumefile_line for each residue in x.
   x = x0 + x1*N + ... + xk*N^k, xi are the residues (this is a hack for GPU)
   FIXME : x0 corresponds to sigma + gpu_curves-1
           xk corresponds to sigma
             should be the other way around

   Returns 1 on success, 0 on error */
int  
write_resumefile (char *fn, int method, mpz_t N, ecm_params params,
		  mpcandi_t *n, mpz_t orig_x0, mpz_t orig_y0, 
		  const char *comment)
{
  FILE *file;
  unsigned int i = 0;
#if defined(HAVE_FCNTL) && defined(HAVE_FILENO)
  struct flock lock;
  int r, fd;
#endif
  mpz_t tmp_x, tmp_y;
  mpz_init (tmp_x);
  mpz_init (tmp_y);

  /* first try to open the file */
#ifdef DEBUG
  if (fn == NULL)
    {
      fprintf (stderr, "write_resumefile: fn == NULL\n");
      exit (EXIT_FAILURE);
    }
#endif
  
  file = fopen (fn, "a");
  if (file == NULL)
    {
      fprintf (stderr, "Could not open file %s for writing\n", fn);
      return 0;
    }
  
#if defined(HAVE_FCNTL) && defined(HAVE_FILENO)
  /* Try to get a lock on the file so several processes can append to
     the same file safely */
  
  /* Supposedly some implementations of fcntl() can get confused over
     garbage in unused fields in a flock struct, so zero it */
  memset (&lock, 0, sizeof (struct flock));
  fd = fileno (file);
  lock.l_type = F_WRLCK;
  lock.l_whence = SEEK_SET;
  lock.l_start = 0;
  lock.l_len = 1; 
  /* F_SETLKW: blocking exclusive lock request */
  r = fcntl (fd, F_SETLKW, &lock);
  if (r != 0)
    {
      fclose (file);
      return 0;
    }

  fseek (file, 0, SEEK_END);
#endif
  

  /* Now can call write_resumefile_line to write in the file */
  if (params->gpu == 0)
    {
      /* Reduce stage 1 residue wrt new cofactor, in case a factor was 
         found */
      mpz_mod (tmp_x, params->x, n->n); 

      /* We write the B1done value to the save file. This requires that
         a correct B1done is returned by the factoring functions. */
      /* FIXME: clang says that params->y == NULL is always false. */
      if (params->y == NULL)
	{
	  write_resumefile_line (file, method, params->B1done, params->sigma,
				 params->sigma_is_A, params->E->type, 
				 params->param, 
				 tmp_x, NULL, n, orig_x0, orig_y0,
				 comment);
	}
      else
	{
	  mpz_mod (tmp_y, params->y, n->n);
	  write_resumefile_line (file, method, params->B1done, params->sigma,
				 params->sigma_is_A, params->E->type,
				 params->param, 
				 tmp_x, tmp_y, n, orig_x0, orig_y0,
				 comment);
	}
    }
  else
    {
      mpz_add_ui (params->sigma, params->sigma, params->gpu_number_of_curves);
      for (i = 0; i < params->gpu_number_of_curves; i++)
        {
          mpz_sub_ui (params->sigma, params->sigma, 1);
          mpz_fdiv_qr (params->x, tmp_x, params->x, N); 
          mpz_mod (tmp_x, tmp_x, n->n);
          write_resumefile_line (file, method, params->B1done, params->sigma,
				 params->sigma_is_A, params->E->type,
				 params->param, 
				 tmp_x, NULL, n, orig_x0, orig_y0, 
				 comment);
        }
    }

  /* closing the file */
#if defined(HAVE_FCNTL) && defined(HAVE_FILENO)
  lock.l_type = F_UNLCK;
  lock.l_whence = SEEK_SET;
  lock.l_start = 0;
  lock.l_len = 1;  
  fcntl (fd, F_SETLKW, &lock); /* F_SETLKW: blocking lock request */
#endif
  fclose (file);

  mpz_clear (tmp_x);
  mpz_clear (tmp_y);

  return 0;
}


/* For the batch mode */
/* Write the batch exponent s in a file */
/* Return the number of bytes written */
int
write_s_in_file (char *fn, mpz_t s)
{
  FILE *file;
  int ret = 0;

#ifdef DEBUG
  if (fn == NULL)
    {
      fprintf (stderr, "write_s_in_file: fn == NULL\n");
      exit (EXIT_FAILURE);
    }
#endif
  
  file = fopen (fn, "wb");
  if (file == NULL)
    {
      fprintf (stderr, "Could not open file %s for writing\n", fn);
      return 0;
    }
  
  ret = mpz_out_raw (file, s);
  
  fclose (file);
  return ret;
}

/* For the batch mode */
/* read the batch exponent s from a file */
int
read_s_from_file (mpz_t s, char *fn, double B1) 
{
  FILE *file;
  mpz_t tmp, tmp2;
  unsigned int val2;
  int ret = 0;

#ifdef DEBUG
  if (fn == NULL)
    {
      fprintf (stderr, "read_s_from_file: fn == NULL\n");
      return 1;
    }
#endif
  
  file = fopen (fn, "rb");
  if (file == NULL)
    {
      fprintf (stderr, "Could not open file %s for reading\n", fn);
      return 1;
    }
 
  ret = mpz_inp_raw (s, file);
  if (ret == 0)
    {
      fprintf (stderr, "read_s_from_file: 0 bytes read from %s\n", fn);
      return 1;
    }

  fclose (file);
          
  /* Some elementaty check that it correspond to the actual B1 */
  mpz_init (tmp);
  mpz_init (tmp2);
  /* check that the valuation of 2 is correct */
  val2 = mpz_scan1 (s, 0);
  mpz_ui_pow_ui (tmp, 2, val2);
  mpz_ui_pow_ui (tmp2, 2, val2+1);
  if (mpz_cmp_d (tmp, B1) > 0 || mpz_cmp_d (tmp2, B1) <= 0)
    {
      fprintf (stderr, "Error, the value of the batch product in %s "
               "does not correspond to B1=%1.0f.\n", fn, B1);
      return 1;
    }

  /* Check that next_prime (B1) does not divide batch_s */
  mpz_set_d (tmp, B1);
  mpz_nextprime (tmp2, tmp);
  if (mpz_divisible_p (s, tmp2))
    {
      fprintf (stderr, "Error, the value of the batch product in %s "
               "does not correspond to B1=%1.0f.\n", fn, B1);
      return 1;
    }

  /* Check that next_prime (sqrt(B1)) divide batch_s only once */
  mpz_set_d (tmp, sqrt(B1));
  mpz_nextprime (tmp2, tmp);
  mpz_mul (tmp, tmp2, tmp2);
  if (!mpz_divisible_p (s, tmp2) || mpz_divisible_p (s, tmp))
    {
      fprintf (stderr, "Error, the value of the batch product in %s "
               "does not correspond to B1=%1.0f.\n", fn, B1);
      return 1;
    }

  mpz_clear (tmp);
  mpz_clear (tmp2);
          
  return 0;
}
