//
//  BenchmarkController.m
//  Prime95
//
//  Created by George Woltman on 1/25/17.
//  Copyright 2017 Mersenne Research, Inc. All rights reserved.
//

#import "BenchmarkController.h"
#include "prime95.h"

@implementation BenchmarkController

- (id)init
{
	if (![super initWithWindowNibName:@"Benchmark"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	[self setBenchType:0];
	[self setMinFFT:IniGetInt (INI_FILE, "MinBenchFFT", 2048)];
	[self setMaxFFT:IniGetInt (INI_FILE, "MaxBenchFFT", 8192)];
	[self setMinBenchableFFT:4];
	[self setMaxBenchableFFT:32768];
	[self setAllFFTSizes:!IniGetInt (INI_FILE, "OnlyBench5678", 1)];
	[self setBenchOneCore:IniGetInt (INI_FILE, "BenchOneCore", 0)];
	[self setBenchAllCores:IniGetInt (INI_FILE, "BenchAllCores", 1)];
	[self setBenchInBetweenCores:IniGetInt (INI_FILE, "BenchInBetweenCores", 0)];
	[self setHyperthreading:IniGetInt (INI_FILE, "BenchHyperthreads", 1)];
	[self setBenchOneWorker:IniGetInt (INI_FILE, "BenchOneWorkerCase", 1)];
	[self setBenchMaxWorkers:IniGetInt (INI_FILE, "BenchMaxWorkersCase", 1)];
	[self setBenchInBetweenWorkers:IniGetInt (INI_FILE, "BenchInBetweenWorkerCases", 0)];
	[self setBenchTime:IniGetInt (INI_FILE, "BenchTime", 15)];
	[self setMinBenchTime:5];
	[self setMaxBenchTime:60];
	[self setEnableds];
}

- (void)setEnableds
{
	[self setMinmaxFFTEnabled: (benchType != 2)];
	[self setAllFFTSizesEnabled: (benchType != 2 && minFFT != maxFFT)];
	[self setBenchOneCoreEnabled: (NUM_CPUS > 1)];
	[self setBenchAllCoresEnabled: (NUM_CPUS > 1)];
	[self setBenchInBetweenCoresEnabled: (NUM_CPUS > 2)];
	[self setHyperthreadingEnabled: (CPU_HYPERTHREADS > 1)];
	[self setBenchOneWorkerEnabled: (benchType == 0 && NUM_CPUS > 1)];
	[self setBenchMaxWorkersEnabled: (benchType == 0 && NUM_CPUS > 1)];
	[self setBenchInBetweenWorkersEnabled: (benchType == 0 && NUM_CPUS > 2)];
	[self setBenchTimeEnabled: (benchType == 0)];
}

- (int)benchType
{
	return benchType;
}

- (void)setBenchType:(int) _value
{
	benchType = _value;
	[self setEnableds];
}

- (int)minFFT
{
	return minFFT;
}

- (void)setMinFFT:(int) _value
{
	minFFT = _value;
	[self setEnableds];
}

- (int)maxFFT
{
	return maxFFT;
}

- (void)setMaxFFT:(int) _value
{
	maxFFT = _value;
	[self setEnableds];
}

@synthesize minmaxFFTEnabled;
@synthesize minBenchableFFT;
@synthesize maxBenchableFFT;
@synthesize allFFTSizes;
@synthesize allFFTSizesEnabled;
@synthesize benchOneCore;
@synthesize benchOneCoreEnabled;
@synthesize benchAllCores;
@synthesize benchAllCoresEnabled;
@synthesize benchInBetweenCores;
@synthesize benchInBetweenCoresEnabled;
@synthesize hyperthreading;
@synthesize hyperthreadingEnabled;
@synthesize benchOneWorker;
@synthesize benchOneWorkerEnabled;
@synthesize benchMaxWorkers;
@synthesize benchMaxWorkersEnabled;
@synthesize benchInBetweenWorkers;
@synthesize benchInBetweenWorkersEnabled;
@synthesize benchTime;
@synthesize benchTimeEnabled;
@synthesize minBenchTime;
@synthesize maxBenchTime;

- (IBAction)ok:(id)sender
{
	[[self window] makeFirstResponder:nil];			// End any active text field edits

	IniWriteInt (INI_FILE, "MinBenchFFT", minFFT);
	IniWriteInt (INI_FILE, "MaxBenchFFT", maxFFT);
	IniWriteInt (INI_FILE, "OnlyBench5678", !allFFTSizes);
	IniWriteInt (INI_FILE, "BenchOneCore", benchOneCore);
	IniWriteInt (INI_FILE, "BenchAllCores", benchAllCores);
	IniWriteInt (INI_FILE, "BenchInBetweenCores", benchInBetweenCores);
	IniWriteInt (INI_FILE, "BenchHyperthreads", hyperthreading);
	IniWriteInt (INI_FILE, "BenchOneWorkerCase", benchOneWorker);
	IniWriteInt (INI_FILE, "BenchMaxWorkersCase", benchMaxWorkers);
	IniWriteInt (INI_FILE, "BenchInBetweenWorkerCases", benchInBetweenWorkers);
	IniWriteInt (INI_FILE, "BenchTime", benchTime);
	LaunchBench (benchType);

	[[self window] performClose:self];
}

@end
