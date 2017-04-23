//
//  BenchmarkController.h
//  Prime95
//
//  Created by George Woltman on 1/25/17.
//  Copyright 2017 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface BenchmarkController : NSWindowController {
	int	benchType;
	int	minFFT;
	int	maxFFT;
	int	minmaxFFTEnabled;
	int	minBenchableFFT;
	int	maxBenchableFFT;
	int	allFFTSizes;
	int	allFFTSizesEnabled;
	int	benchOneCore;
	int	benchOneCoreEnabled;
	int	benchAllCores;
	int	benchAllCoresEnabled;
	int	benchInBetweenCores;
	int	benchInBetweenCoresEnabled;
	int	hyperthreading;
	int	hyperthreadingEnabled;
	int	benchOneWorker;
	int	benchOneWorkerEnabled;
	int	benchMaxWorkers;
	int	benchMaxWorkersEnabled;
	int	benchInBetweenWorkers;
	int	benchInBetweenWorkersEnabled;
	int	benchTime;
	int	benchTimeEnabled;
	int	minBenchTime;
	int	maxBenchTime;
}

@property(readwrite, assign) int benchType;
@property(readwrite, assign) int minFFT;
@property(readwrite, assign) int maxFFT;
@property(readwrite, assign) int minmaxFFTEnabled;
@property(readwrite, assign) int minBenchableFFT;
@property(readwrite, assign) int maxBenchableFFT;
@property(readwrite, assign) int allFFTSizes;
@property(readwrite, assign) int allFFTSizesEnabled;
@property(readwrite, assign) int benchOneCore;
@property(readwrite, assign) int benchOneCoreEnabled;
@property(readwrite, assign) int benchAllCores;
@property(readwrite, assign) int benchAllCoresEnabled;
@property(readwrite, assign) int benchInBetweenCores;
@property(readwrite, assign) int benchInBetweenCoresEnabled;
@property(readwrite, assign) int hyperthreading;
@property(readwrite, assign) int hyperthreadingEnabled;
@property(readwrite, assign) int benchOneWorker;
@property(readwrite, assign) int benchOneWorkerEnabled;
@property(readwrite, assign) int benchMaxWorkers;
@property(readwrite, assign) int benchMaxWorkersEnabled;
@property(readwrite, assign) int benchInBetweenWorkers;
@property(readwrite, assign) int benchInBetweenWorkersEnabled;
@property(readwrite, assign) int benchTime;
@property(readwrite, assign) int benchTimeEnabled;
@property(readwrite, assign) int minBenchTime;
@property(readwrite, assign) int maxBenchTime;

	- (void)reInit;
	- (void)setEnableds;
- (IBAction)ok:(id)sender;

@end
