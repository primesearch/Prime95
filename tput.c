#define GWNUM_FFT_IMPL_VERSION	"28.9"		/* This version number changes whenever FFT implementations change. */
						/* By default, gwnum picks a default FFT implementation based on the CPU architecture. */
						/* You can fine-tune performance by benchmarking the different FFT implementations */
						/* on the end user machine looking for an FFT implementation that is faster than the */
						/* default selection.  Any changes to this version number means you must */
						/* discard old benchmarking data for your fine-tuning. */



if (!lldata->gwdata.bench_pick_nth_fft) {
tput_use_larger_FFTlen_if_faster (); ???
tput_use_larger_FFTlen_if_faster_factoring_in_error_checking_cost (); ???
		tput_select_fastest_fft_impl (gwdata, 1.0, 2, p, c);
	}


	// when near an fft crossover we want to compare smaller speed with error checking to the larger fft without error checking!
tput_use_larger_FFTlen_if_faster ();
tput_use_larger_FFTlen_if_faster_factoring_in_error_checking_cost ();
	if (small_fftlen == large_fftlen || large_fftlen == 0) return (0);



	/* We used to name files xxxx.ini.  Unfortunately, Windows backup/restore */
/* thought these files should be restored to their old values when you */
/* rollback a driver and/or some other reasons.  Now we name our files */
/* xxxx.txt.  This routine converts an old ini name to a new txt name. */

void mangleIniFileName (
	const char *inputFileName,
	char	*outputFileName,
	int	*mangled)	/* Return TRUE if output name != input name */
{
	int	len;
	strcpy (outputFileName, inputFileName);
	len = (int) strlen (outputFileName);
	if (len >= 4 &&
	    outputFileName[len-4] == '.' &&
	    (outputFileName[len-3] == 'I' || outputFileName[len-3] == 'i') &&
	    (outputFileName[len-2] == 'N' || outputFileName[len-2] == 'n') &&
	    (outputFileName[len-1] == 'I' || outputFileName[len-1] == 'i')) {
		strcpy (outputFileName+len-3, "txt");
		*mangled = TRUE;
	} else
		*mangled = FALSE;
}


	
/****************************************************************************/
/*       Utility routines to work with benchmark data in local.txt          */
/****************************************************************************/

struct tput_fft_data {
	struct tput_fft_data *next;		/* Next.  Sorted by date */
	struct tput_fft_data *prev;
	int	date;				/* Benchmark date stored as YYYYMMDD integer */
	float	throughput[10];			/* Room for throughput of up to 10 FFT implementations */
};

struct tput_fft_aggregated {
	struct tput_fft_aggregated *next;	/* Next.  Sorted by all-complex, FFTlen, num_threads */
	struct tput_fft_aggregated *prev;
	unsigned long fftlen;			/* FFT length of throughput benchmark */
	int	all_complex;			/* True if this is a benchmark for an all-complex FFT */
	int	num_threads;			/* Number of threads (CPUs) used in the benchmark */
	int	error_check;			/* TRUE if benchmark was run with roundoff error checking enabled */
	struct tput_fft_data *first_data;	/* Ptr to the raw individual throughput benchamrks */
	float	best_throughput;		/* Throughput of best FFT implementation */
	int	best_implementation;		/* The best FFT implementation */
};

static struct tput_fft_aggregated *THROUGHPUT_DATA;	/* Ptr to all the throughput data generated or read from local.txt */

   routines to add one line to tdata		tput_add_entry (fftlen, all-complex, num-threads, date, throughputstring)
   routines to aggregate / find best		after reading it in?  and after new bench data?
   routines to return best fft impl		tput_best_fft_impl (fftlen, all-complex, num_threads)   
   routines to prune tput data		called during timed_event? at read-in?
	   call tput routines		in readinifile, benchmark code, gwsetups
   timed_event to gen more data

/* Read the existing throughput data from local.txt */

void tput_read_ini_file (void)
{
	char	gwnum_version_string[10];
	char	cpuid_brand_string[49];
	char	tput_val[80];
	int	i, j;

/* Get the gwnum version when the benchmark data was created.  If this does not match the current */
/* gwnum version then we must discard the benchmark data (and start regenerating using the current gwnum code). */

	IniSectionGetString (LOCALINI_FILE, "gwnum", "GwnumVersion", gwnum_version_string, sizeof (gwnum_version_string), NULL);
	if (strcmp (gwnum_version_string, GWNUM_FFT_IMPL_VERSION)) return;

/* Get the CPUID brand string when the benchmark data was created.  If this does not match the current */
/* CPUID brand string as may happen when a local.txt is inadvisably copied to a new computer, then we */
/* must discard the benchmark data (and start regenerating using the new CPU). */

	IniSectionGetString (LOCALINI_FILE, "gwnum", "CpuBrand", cpuid_brand_string, sizeof (cpuid_brand_string), NULL);
	if (strcmp (cpuid_brand_string, CPU_BRAND)) return;

/* We can use the existing throughput benchmark data, loop reading it in.  Format for throughput lines is: */
/*	ThroughputData=date,fftsize/all-complex,numthreads,errchk,tput1,tput2... */

	for (i = 1; ; i++) {
		IniSectionGetNthString (LOCALINI_FILE, "gwnum", "ThroughputData", i, tput_val, sizeof (tput_val), NULL);

....		parse string and save it
	}
		
}

/* Write the throughput data to local.txt */

void tput_write_ini_file (void)
{
	char	gwnum_version_string[10];
	char	cpuid_brand_string[49];
	char	tput_val[80];
	int	i, j;
	struct tput_fft_aggregated *tput;
	struct tput_fft_data *tputdata;

/* Write the gwnum version and CPU brand */

	IniSectionWriteString (LOCALINI_FILE, "gwnum", "GwnumVersion", GWNUM_FFT_IMPL_VERSION);
	IniSectionWriteString (LOCALINI_FILE, "gwnum", "CpuBrand", CPU_BRAND);

/* Loop writing out the throughput benchmark data.  But first clear out the existing throughput data.  Format is: */
/*	ThroughputData=date,fftsize/all-complex,numthreads,errchk,tput1,tput2... */

	IniSectionWriteNthString (LOCALINI_FILE, "gwnum", "ThroughputData", 0, NULL);
	for (tput = THROUGHPUT_DATA, i = 1; tput != NULL; tput = tput->next, i++) {
		tputdata = tput->first_data;
		sprintf (tput_val, "%d,%dK%s,%d,%d",
			 tputdata->date, tput->fftlen >> 10, tput->all_complex ? "ac" : "", tput->num_threads, tput->error_check);
		for (j = 0; j < 10 && tput_data->thoughput[j] > 0.0; j++)
			sprintf (tput_val+strlen(tput_val), ",%.2f", tputdata->throughput[j]);
		IniSectionWriteNthString (LOCALINI_FILE, "gwnum", "ThroughputData", i, tput_val);
	}
}

