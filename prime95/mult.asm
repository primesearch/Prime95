; Copyright 1995-2000 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements the setup, common routines, and global variables
; for the various discrete-weighted transforms.
;

	TITLE   setup

	.386

_TEXT32A SEGMENT PAGE USE32 PUBLIC 'DATA'

EXTRN	_CPU_TYPE:DWORD
EXTRN	_PARG:DWORD
EXTRN	_FFTLEN:DWORD
EXTRN	_INFP:DWORD
EXTRN	_INFF:DWORD
EXTRN	_INFT:DWORD
EXTRN	_SRCARG:DWORD
EXTRN	_SRC2ARG:DWORD
EXTRN	_DESTARG:DWORD
EXTRN	_DEST2ARG:DWORD
EXTRN	_GWPROCPTRS:DWORD
EXTRN	_NORM_ARRAY1:DWORD
EXTRN	_NORM_ARRAY2:DWORD
EXTRN	_PLUS1:DWORD
EXTRN	_GWERROR:DWORD

EXTRN	gwadd1:PROC
EXTRN	gwaddq1:PROC
EXTRN	gwsub1:PROC
EXTRN	gwsubq1:PROC
EXTRN	gwaddsub1:PROC
EXTRN	gwaddsubq1:PROC
EXTRN	gwadd2:PROC
EXTRN	gwaddq2:PROC
EXTRN	gwsub2:PROC
EXTRN	gwsubq2:PROC
EXTRN	gwaddsub2:PROC
EXTRN	gwaddsubq2:PROC
EXTRN	gwadd3:PROC
EXTRN	gwaddq3:PROC
EXTRN	gwsub3:PROC
EXTRN	gwsubq3:PROC
EXTRN	gwaddsub3:PROC
EXTRN	gwaddsubq3:PROC
EXTRN	gwadd4:PROC
EXTRN	gwaddq4:PROC
EXTRN	gwsub4:PROC
EXTRN	gwsubq4:PROC
EXTRN	gwaddsub4:PROC
EXTRN	gwaddsubq4:PROC

exfft	MACRO fft_length
	exfft1	fft_length, _1
	exfft1	fft_length, _2
	exfft1	fft_length, _3
	exfft1	fft_length, _4
	ENDM
exfft1	MACRO fft_length, suffix
	EXTRN	fft&fft_length&suffix:PROC
	EXTRN	fft&fft_length&suffix&PPRO:PROC
	ENDM

exfft	32
exfft	40
exfft	48
;exfft	56
exfft	64
exfft	80
exfft	96
;exfft	112
exfft	128
exfft	160
exfft	192
;exfft	224
exfft	256
exfft	320
exfft	384
exfft	448
exfft	512
exfft	640
exfft	768
exfft	896
exfft	1024
exfft	1280
exfft	1536
exfft	1792
exfft	2048
exfft	2560
exfft	3072
exfft	3584
exfft	4096
exfft	5120
exfft	6144
exfft	7168
exfft	8192
exfft	10K
exfft	12K
exfft	14K
exfft	16K
exfft	20K
exfft	24K
exfft	28K
exfft	32K
exfft	40K
exfft	48K
exfft	56K
exfft	64K
exfft	80K
exfft	96K
exfft	112K
exfft	128K
exfft	160K
exfft	192K
exfft	224K
exfft	256K
exfft	320K
exfft	384K
exfft	448K
exfft	512K
exfft	640K
exfft	768K
exfft	896K
exfft	1024K
exfft	1280K
exfft	1536K
exfft	1792K
exfft	2048K
exfft	2560K
exfft	3072K
exfft	3584K
exfft	4096K

exfft	32p
exfft	64p
exfft	128p
exfft	256p
exfft	512p
exfft	1024p
exfft	2048p
exfft	4096p
exfft	8192p
exfft	16Kp
exfft	32Kp
exfft	64Kp
exfft	128Kp
exfft	256Kp
exfft	512Kp
exfft	1024Kp
exfft	2048Kp
exfft	4096Kp

PUBLIC	sincos_real_data
PUBLIC	sincos_complex_data
PUBLIC	SQRTHALF
PUBLIC	HALF
PUBLIC	BIGVAL
PUBLIC	scaled_numbig
PUBLIC	scaled_numlit
PUBLIC	scaling_ff
PUBLIC	scaling_ff2
PUBLIC	ttmp_ff
PUBLIC	ttmp_ff_inv
PUBLIC	P309
PUBLIC	M809
PUBLIC	M262
PUBLIC	M382
PUBLIC	P951
PUBLIC	P588
PUBLIC	M162
PUBLIC	P618
PUBLIC	P623
PUBLIC	M358
PUBLIC	P404
PUBLIC	P975
PUBLIC	P445
PUBLIC	P180
PUBLIC	M223
PUBLIC	M901
PUBLIC	M691
PUBLIC	P866
PUBLIC	P433
PUBLIC	P577
PUBLIC	P25
PUBLIC	P75
PUBLIC	P3
PUBLIC	limit_inverse_high
PUBLIC	limit_inverse_low
PUBLIC	limit_bigmax_high
PUBLIC	limit_bigmax_low
PUBLIC	limit_ttp_mult_high
PUBLIC	limit_ttp_mult_low
PUBLIC	limit_ttmp_mult_high
PUBLIC	limit_ttmp_mult_low
PUBLIC	pass1_premults
PUBLIC	pass2_premults
PUBLIC	plus1_premults
PUBLIC	count1
PUBLIC	count2
PUBLIC	count3
PUBLIC	count4
PUBLIC	sincos1
PUBLIC	sincos2
PUBLIC	sincos3
PUBLIC	sincos4
PUBLIC	extra_bits
PUBLIC	ffttype
PUBLIC	addcount1
PUBLIC	addcount2
PUBLIC	normcount1
PUBLIC	normcount2
PUBLIC	normgrpptr
PUBLIC	addsubtemp

;
; Global variables needed in multiplication routines
;

	align	32
sincos_real_data	DQ 768 DUP (0.0); Sin/cos values used in real FFTs
sincos_complex_data	DQ 1536 DUP (0.0); Sin/cos values used in complex FFTs
SQRTHALF	DQ	0.0		; Used in all ffts
HALF		DD	0.5		; Used in all ffts
BIGVAL		DD	0.0		; Used to round to an integer
scaled_numbig	DD	0		; numbig * (2^32 / n)
scaled_numlit	DD	0		; numlit * (2^32 / n)
scaling_ff	DD	0		; Fudge factor used in normalizing code
scaling_ff2	DD	0		; Fudge factor used in normalizing code
ttmp_ff_inv	DQ	0.0		; Inverse FFT adjust (2/FFTLEN)
ttmp_ff		DD	0.0		; Two-to-minus-phi adjust (FFTLEN/2)
	align 8
P309		DQ	0.309		; Used in five_reals_fft/unfft
M809		DQ	-0.809		; Used in five_reals_fft/unfft
M262		DQ	-2.618		; Used in five_reals_fft/unfft
M382		DQ	-0.382		; Used in five_reals_fft/unfft
P951		DQ	0.951		; Used in five_reals_fft/unfft
P588		DQ	0.588		; Used in five_reals_fft/unfft
M162		DQ	-1.617		; Used in five_reals_fft/unfft
P618		DQ	0.618		; Used in five_reals_fft/unfft
P623		DQ	0.623		; Used in seven_reals_fft/unfft
M358		DQ	-0.358		; Used in seven_reals_fft/unfft
P404		DQ	4.040		; Used in seven_reals_fft/unfft
P975		DQ	0.975		; Used in seven_reals_fft/unfft
P445		DQ	0.445		; Used in seven_reals_fft/unfft
P180		DQ	1.802		; Used in seven_reals_fft/unfft
M223		DQ	-0.223		; Used in seven_reals_fft/unfft
M901		DQ	-0.901		; Used in seven_reals_fft/unfft
M691		DQ	-0.691		; Used in seven_reals_fft/unfft
P866		DQ	0.866		; Used in six_reals_fft/unfft
P433		DQ	0.433		; Used in six_reals_fft/unfft
P577		DQ	0.577		; Used in six_reals_fft/unfft
P25		DD	0.25		; Used in six_reals_fft/unfft
P75		DD	0.75		; Used in six_reals_fft/unfft
P3		DD	3.0		; Used in six_reals_fft/unfft
	align 8
limit_inverse_high	DQ	0.0	; High and low limit inverses
limit_inverse_low	DQ	0.0
limit_inverse		EQU	limit_inverse_low
limit_bigmax_high	DQ	0.0	; High and low limit * BIGVAL - BIGVAL
limit_bigmax_low	DQ	0.0
limit_bigmax		EQU	limit_bigmax_low
limit_ttp_mult_high	DD	1.0	; High and low limit two-to-phi mult
limit_ttp_mult_low	DD	0.5
limit_ttp_mult		EQU	limit_ttp_mult_low
limit_ttmp_mult_high	DD	1.0	; High/low limit two-to-minus-phi mult
limit_ttmp_mult_low	DD	2.0
limit_ttmp_mult		EQU	limit_ttmp_mult_low
plus1_premults	DD	0		; Address of 2^N+1 premultiplier data
pass1_premults	DD	0		; Address of pass 1 premultiplier data
pass2_premults	DD	0		; Address of pass 2 premultiplier data
count1		DD	0		; Counters used in common fft code
count2		DD	0
count3		DD	0
count4		DD	0
sincos1		DD	0		; Pointers to sine/cosine tables
sincos2		DD	0		; used in common fft code
sincos3		DD	0
sincos4		DD	0
extra_bits	DD	0		; Number of unnormalized adds
					; that can be safely performed.
ffttype		DD	0		; Type of fft (1, 2, 3, or 4)
addcount1	DD	0		; Loop counters used in adding
addcount2	DD	0
normcount1	DD	0		; Loop counters used in normalizing
normcount2	DD	0
normgrpptr	DD	0
addsubtemp	DD	0
		; These values only used during setup
numlit		DD	0		; number of small words in the fft
numbig		DD	0		; number of big words in the fft
P5		DD	5.0		; Used in pfa_5_setup
P7		DD	7.0		; Used in pfa_7_setup
NNNN		DD	0		; Used in premultiplier setup
INCR		DD	0		; Used in premultiplier setup
NOVER16		DD	0		; Used in premultiplier setup
NOVER4		DD	0		; Used in premultiplier setup
SZERO1		DD	0		; Used in premultiplier setup
SZERO2		DD	0		; Used in premultiplier setup
GRPS		DD	0		; Used in premultiplier setup
GRPSIZ		DD	0		; Used in premultiplier setup

jmptable DD	753,	32,	0.0000036,	512
	DD			OFFSET fft32_1, OFFSET fft32_2
	DD			OFFSET fft32_3, OFFSET fft32_4
	DD			OFFSET fft32_1PPRO, OFFSET fft32_2PPRO
	DD			OFFSET fft32_3PPRO, OFFSET fft32_4PPRO
	DD			1, 1, 1, 1
	DD			1, 1, 1, 1, 0
	DD	935,	40,	0.0000057,	928
	DD			OFFSET fft40_1, OFFSET fft40_2
	DD			OFFSET fft40_3, OFFSET fft40_4
	DD			OFFSET fft40_1PPRO, OFFSET fft40_2PPRO
	DD			OFFSET fft40_3PPRO, OFFSET fft40_4PPRO
	DD			1, 1, 1, 1
	DD			1, 1, 1, 1, 20, 40, 0
	DD	1107,	48,	0.0000065,	1104
	DD			OFFSET fft48_1, OFFSET fft48_2
	DD			OFFSET fft48_3, OFFSET fft48_4
	DD			OFFSET fft48_1PPRO, OFFSET fft48_2PPRO
	DD			OFFSET fft48_3PPRO, OFFSET fft48_4PPRO
	DD			1, 1, 1, 1
	DD			1, 1, 1, 1, 24, 48, 0
;	DD	1295,	56,	0.0000084,	1328
;	DD			OFFSET fft56_1, OFFSET fft56_2
;	DD			OFFSET fft56_3, OFFSET fft56_4
;	DD			OFFSET fft56_1PPRO, OFFSET fft56_2PPRO
;	DD			OFFSET fft56_3PPRO, OFFSET fft56_4PPRO
;	DD			1, 1, 1, 1
;	DD			1, 1, 1, 1, 28, 56, 0
	DD	1489,	64,	0.0000083,	1024
	DD			OFFSET fft64_1, OFFSET fft64_2
	DD			OFFSET fft64_3, OFFSET fft64_4
	DD			OFFSET fft64_1PPRO, OFFSET fft64_2PPRO
	DD			OFFSET fft64_3PPRO, OFFSET fft64_4PPRO
	DD			1, 1, 1, 1
	DD			1, 1, 1, 1, 0
	DD	1849,	80,	0.0000121,	1808
	DD			OFFSET fft80_1, OFFSET fft80_2
	DD			OFFSET fft80_3, OFFSET fft80_4
	DD			OFFSET fft80_1PPRO, OFFSET fft80_2PPRO
	DD			OFFSET fft80_3PPRO, OFFSET fft80_4PPRO
	DD			1, 1, 1, 1
	DD			2, 8, 1, 1, 20, 80, 0
	DD	2195,	96,	0.0000141,	2160
	DD			OFFSET fft96_1, OFFSET fft96_2
	DD			OFFSET fft96_3, OFFSET fft96_4
	DD			OFFSET fft96_1PPRO, OFFSET fft96_2PPRO
	DD			OFFSET fft96_3PPRO, OFFSET fft96_4PPRO
	DD			1, 1, 1, 1
	DD			2, 10, 1, 1, 24, 96, 0
;	DD	2557,	112,	0.0000179,	2560
;	DD			OFFSET fft112_1, OFFSET fft112_2
;	DD			OFFSET fft112_3, OFFSET fft112_4
;	DD			OFFSET fft112_1PPRO, OFFSET fft112_2PPRO
;	DD			OFFSET fft112_3PPRO, OFFSET fft112_4PPRO
;	DD			1, 1, 1, 1
;	DD			3, 12, 1, 1, 28, 112, 0
	DD	2935,	128,	0.0000178,	2048
	DD			OFFSET fft128_1, OFFSET fft128_2
	DD			OFFSET fft128_3, OFFSET fft128_4
	DD			OFFSET fft128_1PPRO, OFFSET fft128_2PPRO
	DD			OFFSET fft128_3PPRO, OFFSET fft128_4PPRO
	DD			1, 1, 1, 1
	DD			3, 14, 1, 1, 0
	DD	3641,	160,	0.0000296,	3584
	DD			OFFSET fft160_1, OFFSET fft160_2
	DD			OFFSET fft160_3, OFFSET fft160_4
	DD			OFFSET fft160_1PPRO, OFFSET fft160_2PPRO
	DD			OFFSET fft160_3PPRO, OFFSET fft160_4PPRO
	DD			0ffff0000h+40, 1
	DD			0ffff0000h+5*256+4, 1
	DD			2, 8, 16, 1, 20, 80, 160, 0
	DD	4341,	192,	0.000035,	3872
	DD			OFFSET fft192_1, OFFSET fft192_2
	DD			OFFSET fft192_3, OFFSET fft192_4
	DD			OFFSET fft192_1PPRO, OFFSET fft192_2PPRO
	DD			OFFSET fft192_3PPRO, OFFSET fft192_4PPRO
	DD			0ffff0000h+48, 1
	DD			0ffff0000h+6*256+4, 1
	DD			2, 10, 20, 1, 24, 96, 192, 0
;	DD	5049,	224,	0.000045,	4240
;	DD			OFFSET fft224_1, OFFSET fft224_2
;	DD			OFFSET fft224_3, OFFSET fft224_4
;	DD			OFFSET fft224_1PPRO, OFFSET fft224_2PPRO
;	DD			OFFSET fft224_3PPRO, OFFSET fft224_4PPRO
;	DD			0ffff0000h+56, 1
;	DD			0ffff0000h+7*256+4, 1
;	DD			3, 12, 24, 1, 28, 112, 224, 0
	DD	5797,	256,	0.000045,	2176
	DD			OFFSET fft256_1, OFFSET fft256_2
	DD			OFFSET fft256_3, OFFSET fft256_4
	DD			OFFSET fft256_1PPRO, OFFSET fft256_2PPRO
	DD			OFFSET fft256_3PPRO, OFFSET fft256_4PPRO
	DD			0ffff0000h+64, 1
	DD			0ffff0000h+8*256+4, 1
	DD			3, 14, 28, 1, 0
	DD	7211,	320,	0.000062,	4608
	DD			OFFSET fft320_1, OFFSET fft320_2
	DD			OFFSET fft320_3, OFFSET fft320_4
	DD			OFFSET fft320_1PPRO, OFFSET fft320_2PPRO
	DD			OFFSET fft320_3PPRO, OFFSET fft320_4PPRO
	DD			0ffff0000h+80, 1
	DD			0ffff0000h+10*256+4, 1
	DD			2, 9, 39, 1, 20, 80, 320, 0
	DD	8567,	384,	0.000075,	5120
	DD			OFFSET fft384_1, OFFSET fft384_2
	DD			OFFSET fft384_3, OFFSET fft384_4
	DD			OFFSET fft384_1PPRO, OFFSET fft384_2PPRO
	DD			OFFSET fft384_3PPRO, OFFSET fft384_4PPRO
	DD			0ffff0000h+96, 1
	DD			0ffff0000h+12*256+4, 1
	DD			2, 11, 47, 1, 24, 96, 384, 0
	DD	9967,	448,	0.000093,	5680
	DD			OFFSET fft448_1, OFFSET fft448_2
	DD			OFFSET fft448_3, OFFSET fft448_4
	DD			OFFSET fft448_1PPRO, OFFSET fft448_2PPRO
	DD			OFFSET fft448_3PPRO, OFFSET fft448_4PPRO
	DD			0ffff0000h+112, 1
	DD			0ffff0000h+14*256+4, 1
	DD			3, 13, 55, 1, 28, 112, 448, 0
	DD	11469,	512,	0.000097,	2304
	DD			OFFSET fft512_1, OFFSET fft512_2
	DD			OFFSET fft512_3, OFFSET fft512_4
	DD			OFFSET fft512_1PPRO, OFFSET fft512_2PPRO
	DD			OFFSET fft512_3PPRO, OFFSET fft512_4PPRO
	DD			0ffff0000h+128, 1
	DD			0ffff0000h+16*256+4, 1
	DD			3, 15, 63, 1, 0
	DD	14175,	640,	0.000140,	8560
	DD			OFFSET fft640_1, OFFSET fft640_2
	DD			OFFSET fft640_3, OFFSET fft640_4
	DD			OFFSET fft640_1PPRO, OFFSET fft640_2PPRO
	DD			OFFSET fft640_3PPRO, OFFSET fft640_4PPRO
	DD			0ffff0000h+160, 1
	DD			0ffff0000h+20*256+4, 1
	DD			2, 9, 39, 79, 20, 80, 320, 640, 0
	DD	16945,	768,	0.000167,	9872
	DD			OFFSET fft768_1, OFFSET fft768_2
	DD			OFFSET fft768_3, OFFSET fft768_4
	DD			OFFSET fft768_1PPRO, OFFSET fft768_2PPRO
	DD			OFFSET fft768_3PPRO, OFFSET fft768_4PPRO
	DD			0ffff0000h+192, 1
	DD			0ffff0000h+24*256+4, 1
	DD			2, 11, 47, 95, 24, 96, 384, 768, 0
	DD	19739,	896,	0.000210,	11232
	DD			OFFSET fft896_1, OFFSET fft896_2
	DD			OFFSET fft896_3, OFFSET fft896_4
	DD			OFFSET fft896_1PPRO, OFFSET fft896_2PPRO
	DD			OFFSET fft896_3PPRO, OFFSET fft896_4PPRO
	DD			0ffff0000h+224, 1
	DD			0ffff0000h+28*256+4, 1
	DD			3, 13, 55, 111, 28, 112, 448, 896, 0
	DD	22599,	1024,	0.000218,	2560
	DD			OFFSET fft1024_1, OFFSET fft1024_2
	DD			OFFSET fft1024_3, OFFSET fft1024_4
	DD			OFFSET fft1024_1PPRO, OFFSET fft1024_2PPRO
	DD			OFFSET fft1024_3PPRO, OFFSET fft1024_4PPRO
	DD			0ffff0000h+256, 1
	DD			0ffff0000h+32*256+4, 1
	DD			3, 15, 63, 127, 0
	DD	28169,	1280,	0.000302,	6816
	DD			OFFSET fft1280_1, OFFSET fft1280_2
	DD			OFFSET fft1280_3, OFFSET fft1280_4
	DD			OFFSET fft1280_1PPRO, OFFSET fft1280_2PPRO
	DD			OFFSET fft1280_3PPRO, OFFSET fft1280_4PPRO
	DD			1*256+128, 0ffff0000h+64
	DD			1*65536+16*256+4, 0ffff0000h+8*256+4
	DD			9*256+3, 1, 1, 1, 20, 0
	DD	33459,	1536,	0.000365,	7840
	DD			OFFSET fft1536_1, OFFSET fft1536_2
	DD			OFFSET fft1536_3, OFFSET fft1536_4
	DD			OFFSET fft1536_1PPRO, OFFSET fft1536_2PPRO
	DD			OFFSET fft1536_3PPRO, OFFSET fft1536_4PPRO
	DD			1*256+128, 0ffff0000h+128
	DD			1*65536+16*256+4, 0ffff0000h+16*256+4
	DD			11*256+3, 1, 1, 1, 24, 0
	DD	39065,	1792,	0.000456,	8912
	DD			OFFSET fft1792_1, OFFSET fft1792_2
	DD			OFFSET fft1792_3, OFFSET fft1792_4
	DD			OFFSET fft1792_1PPRO, OFFSET fft1792_2PPRO
	DD			OFFSET fft1792_3PPRO, OFFSET fft1792_4PPRO
	DD			2*256+128, 0ffff0000h+64
	DD			2*65536+16*256+4, 0ffff0000h+8*256+4
	DD			13*256+3, 1, 1, 1, 28, 0
	DD	44771,	2048,	0.000490,	9792
	DD			OFFSET fft2048_1, OFFSET fft2048_2
	DD			OFFSET fft2048_3, OFFSET fft2048_4
	DD			OFFSET fft2048_1PPRO, OFFSET fft2048_2PPRO
	DD			OFFSET fft2048_3PPRO, OFFSET fft2048_4PPRO
	DD			2*256+128, 0ffff0000h+128
	DD			2*65536+16*256+4, 0ffff0000h+16*256+4
	DD			15*256+3, 1, 1, 1, 0
	DD	55541,	2560,	0.000708,	12128
	DD			OFFSET fft2560_1, OFFSET fft2560_2
	DD			OFFSET fft2560_3, OFFSET fft2560_4
	DD			OFFSET fft2560_1PPRO, OFFSET fft2560_2PPRO
	DD			OFFSET fft2560_3PPRO, OFFSET fft2560_4PPRO
	DD			3*256+128, 0ffff0000h+128
	DD			3*65536+16*256+4, 0ffff0000h+16*256+4
	DD			19*256+3, 1, 1, 1, 20, 40, 0
	DD	66161,	3072,	0.000851,	14224
	DD			OFFSET fft3072_1, OFFSET fft3072_2
	DD			OFFSET fft3072_3, OFFSET fft3072_4
	DD			OFFSET fft3072_1PPRO, OFFSET fft3072_2PPRO
	DD			OFFSET fft3072_3PPRO, OFFSET fft3072_4PPRO
	DD			4*256+128, 0ffff0000h+128
	DD			4*65536+16*256+4, 0ffff0000h+16*256+4
	DD			23*256+3, 1, 1, 1, 24, 48, 0
	DD	77000,	3584,	0.00107,	16368
	DD			OFFSET fft3584_1, OFFSET fft3584_2
	DD			OFFSET fft3584_3, OFFSET fft3584_4
	DD			OFFSET fft3584_1PPRO, OFFSET fft3584_2PPRO
	DD			OFFSET fft3584_3PPRO, OFFSET fft3584_4PPRO
	DD			5*256+128, 0ffff0000h+128
	DD			5*65536+16*256+4, 0ffff0000h+16*256+4
	DD			27*256+3, 1, 1, 1, 28, 56, 0
	DD	88500,	4096,	0.00113,	17984
	DD			OFFSET fft4096_1, OFFSET fft4096_2
	DD			OFFSET fft4096_3, OFFSET fft4096_4
	DD			OFFSET fft4096_1PPRO, OFFSET fft4096_2PPRO
	DD			OFFSET fft4096_3PPRO, OFFSET fft4096_4PPRO
	DD			6*256+128, 0ffff0000h+128
	DD			6*65536+16*256+4, 0ffff0000h+16*256+4
	DD			31*256+3, 1, 1, 1, 0
	DD	110000,	5120,	0.00152,	22608
	DD			OFFSET fft5120_1, OFFSET fft5120_2
	DD			OFFSET fft5120_3, OFFSET fft5120_4
	DD			OFFSET fft5120_1PPRO, OFFSET fft5120_2PPRO
	DD			OFFSET fft5120_3PPRO, OFFSET fft5120_4PPRO
	DD			8*256+128, 0ffff0000h+128
	DD			8*65536+16*256+4, 0ffff0000h+16*256+4
	DD			39*256+3, 1, 1, 1, 20, 80, 0
	DD	130400,	6144,	0.00191,	26800
	DD			OFFSET fft6144_1, OFFSET fft6144_2
	DD			OFFSET fft6144_3, OFFSET fft6144_4
	DD			OFFSET fft6144_1PPRO, OFFSET fft6144_2PPRO
	DD			OFFSET fft6144_3PPRO, OFFSET fft6144_4PPRO
	DD			10*256+128, 0ffff0000h+128
	DD			10*65536+16*256+4, 0ffff0000h+16*256+4
	DD			47*256+3, 1, 1, 1, 24, 96, 0
	DD	152200,	7168,	0.00226,	31040
	DD			OFFSET fft7168_1, OFFSET fft7168_2
	DD			OFFSET fft7168_3, OFFSET fft7168_4
	DD			OFFSET fft7168_1PPRO, OFFSET fft7168_2PPRO
	DD			OFFSET fft7168_3PPRO, OFFSET fft7168_4PPRO
	DD			12*256+128, 0ffff0000h+128
	DD			12*65536+16*256+4, 0ffff0000h+16*256+4
	DD			55*256+3, 1, 1, 1, 28, 112, 0
	DD	174600,	8192,	0.00242,	34368
	DD			OFFSET fft8192_1, OFFSET fft8192_2
	DD			OFFSET fft8192_3, OFFSET fft8192_4
	DD			OFFSET fft8192_1PPRO, OFFSET fft8192_2PPRO
	DD			OFFSET fft8192_3PPRO, OFFSET fft8192_4PPRO
	DD			14*256+128, 0ffff0000h+128
	DD			14*65536+16*256+4, 0ffff0000h+16*256+4
	DD			63*256+3, 1, 1, 1, 0
	DD	217000,	10240,	0.00333,	26720
	DD			OFFSET fft10K_1, OFFSET fft10K_2
	DD			OFFSET fft10K_3, OFFSET fft10K_4
	DD			OFFSET fft10K_1PPRO, OFFSET fft10K_2PPRO
	DD			OFFSET fft10K_3PPRO, OFFSET fft10K_4PPRO
	DD			19*65536+8*256+16, 1
	DD			19*65536+8*256+4, 1
	DD			19*256+15, 1, 1, 1, 20, 40, 0
	DD	258400,	12288,	0.00397,	31888
	DD			OFFSET fft12K_1, OFFSET fft12K_2
	DD			OFFSET fft12K_3, OFFSET fft12K_4
	DD			OFFSET fft12K_1PPRO, OFFSET fft12K_2PPRO
	DD			OFFSET fft12K_3PPRO, OFFSET fft12K_4PPRO
	DD			23*65536+8*256+16, 1
	DD			23*65536+8*256+4, 1
	DD			23*256+15, 1, 1, 1, 24, 48, 0
	DD	301000,	14336,	0.00488,	37104
	DD			OFFSET fft14K_1, OFFSET fft14K_2
	DD			OFFSET fft14K_3, OFFSET fft14K_4
	DD			OFFSET fft14K_1PPRO, OFFSET fft14K_2PPRO
	DD			OFFSET fft14K_3PPRO, OFFSET fft14K_4PPRO
	DD			27*65536+8*256+16, 1
	DD			27*65536+8*256+4, 1
	DD			27*256+15, 1, 1, 1, 28, 56, 0
	DD	345400,	16384,	0.00522,	41792
	DD			OFFSET fft16K_1, OFFSET fft16K_2
	DD			OFFSET fft16K_3, OFFSET fft16K_4
	DD			OFFSET fft16K_1PPRO, OFFSET fft16K_2PPRO
	DD			OFFSET fft16K_3PPRO, OFFSET fft16K_4PPRO
	DD			31*65536+8*256+16, 1
	DD			31*65536+8*256+4, 1
	DD			31*256+15, 1, 1, 1, 0
	DD	428200,	20480,	0.00692,	52560
	DD			OFFSET fft20K_1, OFFSET fft20K_2
	DD			OFFSET fft20K_3, OFFSET fft20K_4
	DD			OFFSET fft20K_1PPRO, OFFSET fft20K_2PPRO
	DD			OFFSET fft20K_3PPRO, OFFSET fft20K_4PPRO
	DD			39*65536+8*256+16, 1
	DD			39*65536+8*256+4, 1
	DD			39*256+15, 1, 1, 1, 20, 80, 0
	DD	509200,	24576,	0.00826,	62896
	DD			OFFSET fft24K_1, OFFSET fft24K_2
	DD			OFFSET fft24K_3, OFFSET fft24K_4
	DD			OFFSET fft24K_1PPRO, OFFSET fft24K_2PPRO
	DD			OFFSET fft24K_3PPRO, OFFSET fft24K_4PPRO
	DD			47*65536+8*256+16, 1
	DD			47*65536+8*256+4, 1
	DD			47*256+15, 1, 1, 1, 24, 96, 0
	DD	594200,	28672,	0.0101,		73280
	DD			OFFSET fft28K_1, OFFSET fft28K_2
	DD			OFFSET fft28K_3, OFFSET fft28K_4
	DD			OFFSET fft28K_1PPRO, OFFSET fft28K_2PPRO
	DD			OFFSET fft28K_3PPRO, OFFSET fft28K_4PPRO
	DD			55*65536+8*256+16, 1
	DD			55*65536+8*256+4, 1
	DD			55*256+15, 1, 1, 1, 28, 112, 0
	DD	680000,	32768,	0.0109,		82752
	DD			OFFSET fft32K_1, OFFSET fft32K_2
	DD			OFFSET fft32K_3, OFFSET fft32K_4
	DD			OFFSET fft32K_1PPRO, OFFSET fft32K_2PPRO
	DD			OFFSET fft32K_3PPRO, OFFSET fft32K_4PPRO
	DD			63*65536+8*256+16, 1
	DD			63*65536+8*256+4, 1
	DD			63*256+15, 1, 1, 1, 0
	DD	843600,	40960,	0.0151,		104672
	DD			OFFSET fft40K_1, OFFSET fft40K_2
	DD			OFFSET fft40K_3, OFFSET fft40K_4
	DD			OFFSET fft40K_1PPRO, OFFSET fft40K_2PPRO
	DD			OFFSET fft40K_3PPRO, OFFSET fft40K_4PPRO
	DD			79*65536+8*256+16, 1
	DD			79*65536+8*256+4, 1
	DD			79*256+15, 1, 1, 1, 20, 80, 160, 0
	DD	1006400, 49152,	0.0184	,	125440
	DD			OFFSET fft48K_1, OFFSET fft48K_2
	DD			OFFSET fft48K_3, OFFSET fft48K_4
	DD			OFFSET fft48K_1PPRO, OFFSET fft48K_2PPRO
	DD			OFFSET fft48K_3PPRO, OFFSET fft48K_4PPRO
	DD			95*65536+8*256+16, 1
	DD			95*65536+8*256+4, 1
	DD			95*256+15, 1, 1, 1, 24, 96, 192, 0
	DD	1176200, 57344,	0.0227,		146256
	DD			OFFSET fft56K_1, OFFSET fft56K_2
	DD			OFFSET fft56K_3, OFFSET fft56K_4
	DD			OFFSET fft56K_1PPRO, OFFSET fft56K_2PPRO
	DD			OFFSET fft56K_3PPRO, OFFSET fft56K_4PPRO
	DD			111*65536+8*256+16, 1
	DD			111*65536+8*256+4, 1
	DD			111*256+15, 1, 1, 1, 28, 112, 224, 0
	DD	1345000, 65536,	0.0252,		164672
	DD			OFFSET fft64K_1, OFFSET fft64K_2
	DD			OFFSET fft64K_3, OFFSET fft64K_4
	DD			OFFSET fft64K_1PPRO, OFFSET fft64K_2PPRO
	DD			OFFSET fft64K_3PPRO, OFFSET fft64K_4PPRO
	DD			127*65536+8*256+16, 1
	DD			127*65536+8*256+4, 1
	DD			127*256+15, 1, 1, 1, 0
	DD	1665000, 81920, 0.0360,		206528
	DD			OFFSET fft80K_1, OFFSET fft80K_2
	DD			OFFSET fft80K_3, OFFSET fft80K_4
	DD			OFFSET fft80K_1PPRO, OFFSET fft80K_2PPRO
	DD			OFFSET fft80K_3PPRO, OFFSET fft80K_4PPRO
	DD			1, 1
	DD			159*32768, 1
	DD			4*65536+8*256+256/2+16, 1, 1, 1, 0
	DD	1990000, 98304, 0.0445,		247488
	DD			OFFSET fft96K_1, OFFSET fft96K_2
	DD			OFFSET fft96K_3, OFFSET fft96K_4
	DD			OFFSET fft96K_1PPRO, OFFSET fft96K_2PPRO
	DD			OFFSET fft96K_3PPRO, OFFSET fft96K_4PPRO
	DD			1, 1
	DD			191*32768, 1
	DD			5*65536+8*256+256/2+16, 1, 1, 1, 0
	DD	2323000, 114688, 0.0548,	288896
	DD			OFFSET fft112K_1, OFFSET fft112K_2
	DD			OFFSET fft112K_3, OFFSET fft112K_4
	DD			OFFSET fft112K_1PPRO, OFFSET fft112K_2PPRO
	DD			OFFSET fft112K_3PPRO, OFFSET fft112K_4PPRO
	DD			1, 1
	DD			223*32768, 1
	DD			6*65536+8*256+256/2+16, 1, 1, 1, 0
	DD	2655500, 131072, 0.0604,	329856
	DD			OFFSET fft128K_1, OFFSET fft128K_2
	DD			OFFSET fft128K_3, OFFSET fft128K_4
	DD			OFFSET fft128K_1PPRO, OFFSET fft128K_2PPRO
	DD			OFFSET fft128K_3PPRO, OFFSET fft128K_4PPRO
	DD			1, 1
	DD			255*32768, 1
	DD			7*65536+8*256+256/2+16, 1, 1, 1, 0
	DD	3290000, 163840, 0.0830,	411840
	DD			OFFSET fft160K_1, OFFSET fft160K_2
	DD			OFFSET fft160K_3, OFFSET fft160K_4
	DD			OFFSET fft160K_1PPRO, OFFSET fft160K_2PPRO
	DD			OFFSET fft160K_3PPRO, OFFSET fft160K_4PPRO
	DD			1, 1
	DD			319*32768, 1
	DD			9*65536+8*256+256/2+16, 2, 1, 1, 0
	DD	3935000, 196608, 0.0982,	493760
	DD			OFFSET fft192K_1, OFFSET fft192K_2
	DD			OFFSET fft192K_3, OFFSET fft192K_4
	DD			OFFSET fft192K_1PPRO, OFFSET fft192K_2PPRO
	DD			OFFSET fft192K_3PPRO, OFFSET fft192K_4PPRO
	DD			1, 1
	DD			383*32768, 1
	DD			11*65536+8*256+256/2+16, 2, 1, 1, 0
	DD	4598000, 229376, 0.1193,	576384
	DD			OFFSET fft224K_1, OFFSET fft224K_2
	DD			OFFSET fft224K_3, OFFSET fft224K_4
	DD			OFFSET fft224K_1PPRO, OFFSET fft224K_2PPRO
	DD			OFFSET fft224K_3PPRO, OFFSET fft224K_4PPRO
	DD			1, 1
	DD			447*32768, 1
	DD			13*65536+8*256+256/2+16, 2, 1, 1, 0
	DD	5250000, 262144, 0.1316,	658304
	DD			OFFSET fft256K_1, OFFSET fft256K_2
	DD			OFFSET fft256K_3, OFFSET fft256K_4
	DD			OFFSET fft256K_1PPRO, OFFSET fft256K_2PPRO
	DD			OFFSET fft256K_3PPRO, OFFSET fft256K_4PPRO
	DD			1, 1
	DD			511*32768, 1
	DD			15*65536+8*256+256/2+16, 2, 1, 1, 0
	DD	6515000, 327680, 0.1726,	824160
	DD			OFFSET fft320K_1, OFFSET fft320K_2
	DD			OFFSET fft320K_3, OFFSET fft320K_4
	DD			OFFSET fft320K_1PPRO, OFFSET fft320K_2PPRO
	DD			OFFSET fft320K_3PPRO, OFFSET fft320K_4PPRO
	DD			1, 1
	DD			639*32768, 1
	DD			19*65536+8*256+256/2+16, 1, 1, 1, 20, 0
	DD	7730000, 393216, 0.2107,	988896
	DD			OFFSET fft384K_1, OFFSET fft384K_2
	DD			OFFSET fft384K_3, OFFSET fft384K_4
	DD			OFFSET fft384K_1PPRO, OFFSET fft384K_2PPRO
	DD			OFFSET fft384K_3PPRO, OFFSET fft384K_4PPRO
	DD			1, 1
	DD			767*32768, 1
	DD			23*65536+8*256+256/2+16, 1, 1, 1, 24, 0
	DD	9020000, 458752, 0.2520,	1153680
	DD			OFFSET fft448K_1, OFFSET fft448K_2
	DD			OFFSET fft448K_3, OFFSET fft448K_4
	DD			OFFSET fft448K_1PPRO, OFFSET fft448K_2PPRO
	DD			OFFSET fft448K_3PPRO, OFFSET fft448K_4PPRO
	DD			1, 1
	DD			895*32768, 1
	DD			27*65536+8*256+256/2+16, 1, 1, 1, 28, 0
	DD	10320000, 524288, 0.2808,	1318272
	DD			OFFSET fft512K_1, OFFSET fft512K_2
	DD			OFFSET fft512K_3, OFFSET fft512K_4
	DD			OFFSET fft512K_1PPRO, OFFSET fft512K_2PPRO
	DD			OFFSET fft512K_3PPRO, OFFSET fft512K_4PPRO
	DD			1, 1
	DD			1023*32768, 1
	DD			31*65536+8*256+256/2+16, 1, 1, 1, 0
	DD	12830000, 655360, 0.372,	1645664
	DD			OFFSET fft640K_1, OFFSET fft640K_2
	DD			OFFSET fft640K_3, OFFSET fft640K_4
	DD			OFFSET fft640K_1PPRO, OFFSET fft640K_2PPRO
	DD			OFFSET fft640K_3PPRO, OFFSET fft640K_4PPRO
	DD			1, 1
	DD			1279*32768, 1
	DD			39*65536+8*256+256/2+16, 2, 1, 1, 20, 0
	DD	15270000, 786432, 0.453,	1974752
	DD			OFFSET fft768K_1, OFFSET fft768K_2
	DD			OFFSET fft768K_3, OFFSET fft768K_4
	DD			OFFSET fft768K_1PPRO, OFFSET fft768K_2PPRO
	DD			OFFSET fft768K_3PPRO, OFFSET fft768K_4PPRO
	DD			1, 1
	DD			1535*32768, 1
	DD			47*65536+8*256+256/2+16, 2, 1, 1, 24, 0
	DD	17850000, 917504, 0.536,	2303888
	DD			OFFSET fft896K_1, OFFSET fft896K_2
	DD			OFFSET fft896K_3, OFFSET fft896K_4
	DD			OFFSET fft896K_1PPRO, OFFSET fft896K_2PPRO
	DD			OFFSET fft896K_3PPRO, OFFSET fft896K_4PPRO
	DD			1, 1
	DD			1791*32768, 1
	DD			55*65536+8*256+256/2+16, 2, 1, 1, 28, 0
	DD	20400000, 1048576, 0.600,	2632832
	DD			OFFSET fft1024K_1, OFFSET fft1024K_2
	DD			OFFSET fft1024K_3, OFFSET fft1024K_4
	DD			OFFSET fft1024K_1PPRO, OFFSET fft1024K_2PPRO
	DD			OFFSET fft1024K_3PPRO, OFFSET fft1024K_4PPRO
	DD			1, 1
	DD			2047*32768, 1
	DD			63*65536+8*256+256/2+16, 2, 1, 1, 0
	DD	25330000, 1310720, 0.776,	3295632
	DD			OFFSET fft1280K_1, OFFSET fft1280K_2
	DD			OFFSET fft1280K_3, OFFSET fft1280K_4
	DD			OFFSET fft1280K_1PPRO, OFFSET fft1280K_2PPRO
	DD			OFFSET fft1280K_3PPRO, OFFSET fft1280K_4PPRO
	DD			1, 1
	DD			2559*32768, 1
	DD			79*65536+8*256+256/2+16, 1, 1, 1
	DD			20, 80, 0
	DD	30100000, 1572864, 0.934,	3954672
	DD			OFFSET fft1536K_1, OFFSET fft1536K_2
	DD			OFFSET fft1536K_3, OFFSET fft1536K_4
	DD			OFFSET fft1536K_1PPRO, OFFSET fft1536K_2PPRO
	DD			OFFSET fft1536K_3PPRO, OFFSET fft1536K_4PPRO
	DD			1, 1
	DD			3071*32768, 1
	DD			95*65536+8*256+256/2+16, 1, 1, 1
	DD			24, 96, 0
	DD	35100000, 1835008, 1.113,	4613760
	DD			OFFSET fft1792K_1, OFFSET fft1792K_2
	DD			OFFSET fft1792K_3, OFFSET fft1792K_4
	DD			OFFSET fft1792K_1PPRO, OFFSET fft1792K_2PPRO
	DD			OFFSET fft1792K_3PPRO, OFFSET fft1792K_4PPRO
	DD			1, 1
	DD			3583*32768, 1
	DD			111*65536+8*256+256/2+16, 1, 1, 1
	DD			28, 112, 0
	DD	40250000, 2097152, 1.226,	5271936
	DD			OFFSET fft2048K_1, OFFSET fft2048K_2
	DD			OFFSET fft2048K_3, OFFSET fft2048K_4
	DD			OFFSET fft2048K_1PPRO, OFFSET fft2048K_2PPRO
	DD			OFFSET fft2048K_3PPRO, OFFSET fft2048K_4PPRO
	DD			1, 1
	DD			4095*32768, 1
	DD			127*65536+8*256+256/2+16, 1, 1, 1, 0
	DD	49900000, 2621440, 1.636,	6582416
	DD			OFFSET fft2560K_1, OFFSET fft2560K_2
	DD			OFFSET fft2560K_3, OFFSET fft2560K_4
	DD			OFFSET fft2560K_1PPRO, OFFSET fft2560K_2PPRO
	DD			OFFSET fft2560K_3PPRO, OFFSET fft2560K_4PPRO
	DD			1, 1
	DD			5119*32768, 1
	DD			159*65536+8*256+256/2+16, 2, 1, 1
	DD			20, 80, 0
	DD	59400000, 3145728, 1.990,	7898864
	DD			OFFSET fft3072K_1, OFFSET fft3072K_2
	DD			OFFSET fft3072K_3, OFFSET fft3072K_4
	DD			OFFSET fft3072K_1PPRO, OFFSET fft3072K_2PPRO
	DD			OFFSET fft3072K_3PPRO, OFFSET fft3072K_4PPRO
	DD			1, 1
	DD			6143*32768, 1
	DD			191*65536+8*256+256/2+16, 2, 1, 1
	DD			24, 96, 0
	DD	69000000, 3670016, 2.380,	9215360
	DD			OFFSET fft3584K_1, OFFSET fft3584K_2
	DD			OFFSET fft3584K_3, OFFSET fft3584K_4
	DD			OFFSET fft3584K_1PPRO, OFFSET fft3584K_2PPRO
	DD			OFFSET fft3584K_3PPRO, OFFSET fft3584K_4PPRO
	DD			1, 1
	DD			7167*32768, 1
	DD			223*65536+8*256+256/2+16, 2, 1, 1
	DD			28, 112, 0
	DD	79300000, 4194304, 2.604,	10530944
	DD			OFFSET fft4096K_1, OFFSET fft4096K_2
	DD			OFFSET fft4096K_3, OFFSET fft4096K_4
	DD			OFFSET fft4096K_1PPRO, OFFSET fft4096K_2PPRO
	DD			OFFSET fft4096K_3PPRO, OFFSET fft4096K_4PPRO
	DD			1, 1
	DD			8191*32768, 1
	DD			255*65536+8*256+256/2+16, 2, 1, 1, 0
jmptablep DD	753,	32,	0.000004,	768
	DD			OFFSET fft32p_1, OFFSET fft32p_2
	DD			OFFSET fft32p_3, OFFSET fft32p_4
	DD			OFFSET fft32p_1PPRO, OFFSET fft32p_2PPRO
	DD			OFFSET fft32p_3PPRO, OFFSET fft32p_4PPRO
	DD			1, 1, 1, 1
	DD			1, 1, 1, 1, 0
	DD	1489,	64,	0.000010,	1536
	DD			OFFSET fft64p_1, OFFSET fft64p_2
	DD			OFFSET fft64p_3, OFFSET fft64p_4
	DD			OFFSET fft64p_1PPRO, OFFSET fft64p_2PPRO
	DD			OFFSET fft64p_3PPRO, OFFSET fft64p_4PPRO
	DD			1, 1, 1, 1
	DD			1, 1, 1, 1, 0
	DD	2935,	128,	0.000021,	3072
	DD			OFFSET fft128p_1, OFFSET fft128p_2
	DD			OFFSET fft128p_3, OFFSET fft128p_4
	DD			OFFSET fft128p_1PPRO, OFFSET fft128p_2PPRO
	DD			OFFSET fft128p_3PPRO, OFFSET fft128p_4PPRO
	DD			1, 1, 1, 1
	DD			3, 14, 1, 1, 0
	DD	5797,	256,	0.000051,	4224
	DD			OFFSET fft256p_1, OFFSET fft256p_2
	DD			OFFSET fft256p_3, OFFSET fft256p_4
	DD			OFFSET fft256p_1PPRO, OFFSET fft256p_2PPRO
	DD			OFFSET fft256p_3PPRO, OFFSET fft256p_4PPRO
	DD			0ffff0000h+64, 1
	DD			0ffff0000h+8*256+4, 1
	DD			3, 14, 28, 1, 0
	DD	11469,	512,	0.000106,	6400
	DD			OFFSET fft512p_1, OFFSET fft512p_2
	DD			OFFSET fft512p_3, OFFSET fft512p_4
	DD			OFFSET fft512p_1PPRO, OFFSET fft512p_2PPRO
	DD			OFFSET fft512p_3PPRO, OFFSET fft512p_4PPRO
	DD			0ffff0000h+128, 1
	DD			0ffff0000h+16*256+4, 1
	DD			3, 15, 63, 1, 0
	DD	22599,	1024,	0.000249,	10752
	DD			OFFSET fft1024p_1, OFFSET fft1024p_2
	DD			OFFSET fft1024p_3, OFFSET fft1024p_4
	DD			OFFSET fft1024p_1PPRO, OFFSET fft1024p_2PPRO
	DD			OFFSET fft1024p_3PPRO, OFFSET fft1024p_4PPRO
	DD			0ffff0000h+256, 1
	DD			0ffff0000h+32*256+4, 1
	DD			3, 15, 63, 127, 0
	DD	44771,	2048,	0.000582,	26624
	DD			OFFSET fft2048p_1, OFFSET fft2048p_2
	DD			OFFSET fft2048p_3, OFFSET fft2048p_4
	DD			OFFSET fft2048p_1PPRO, OFFSET fft2048p_2PPRO
	DD			OFFSET fft2048p_3PPRO, OFFSET fft2048p_4PPRO
	DD			2*256+128, 0ffff0000h+128
	DD			2*65536+16*256+4, 0ffff0000h+16*256+4
	DD			16*256+4, 1, 1, 1, 0
	DD	88500,	4096,	0.00135,	51200
	DD			OFFSET fft4096p_1, OFFSET fft4096p_2
	DD			OFFSET fft4096p_3, OFFSET fft4096p_4
	DD			OFFSET fft4096p_1PPRO, OFFSET fft4096p_2PPRO
	DD			OFFSET fft4096p_3PPRO, OFFSET fft4096p_4PPRO
	DD			6*256+128, 0ffff0000h+128
	DD			6*65536+16*256+4, 0ffff0000h+16*256+4
	DD			32*256+4, 1, 1, 1, 0
	DD	174600,	8192,	0.00284,	100352
	DD			OFFSET fft8192p_1, OFFSET fft8192p_2
	DD			OFFSET fft8192p_3, OFFSET fft8192p_4
	DD			OFFSET fft8192p_1PPRO, OFFSET fft8192p_2PPRO
	DD			OFFSET fft8192p_3PPRO, OFFSET fft8192p_4PPRO
	DD			14*256+128, 0ffff0000h+128
	DD			14*65536+16*256+4, 0ffff0000h+16*256+4
	DD			64*256+4, 1, 1, 1, 0
	DD	345400,	16384,	0.00588,	174080
	DD			OFFSET fft16Kp_1, OFFSET fft16Kp_2
	DD			OFFSET fft16Kp_3, OFFSET fft16Kp_4
	DD			OFFSET fft16Kp_1PPRO, OFFSET fft16Kp_2PPRO
	DD			OFFSET fft16Kp_3PPRO, OFFSET fft16Kp_4PPRO
	DD			31*65536+8*256+16, 1
	DD			31*65536+8*256+4, 1
	DD			32*256+16, 1, 1, 1, 0
	DD	680000,	32768,	0.01299,	346112
	DD			OFFSET fft32Kp_1, OFFSET fft32Kp_2
	DD			OFFSET fft32Kp_3, OFFSET fft32Kp_4
	DD			OFFSET fft32Kp_1PPRO, OFFSET fft32Kp_2PPRO
	DD			OFFSET fft32Kp_3PPRO, OFFSET fft32Kp_4PPRO
	DD			63*65536+8*256+16, 1
	DD			63*65536+8*256+4, 1
	DD			64*256+16, 1, 1, 1, 0
	DD	1345000, 65536,	0.03283,	690176
	DD			OFFSET fft64Kp_1, OFFSET fft64Kp_2
	DD			OFFSET fft64Kp_3, OFFSET fft64Kp_4
	DD			OFFSET fft64Kp_1PPRO, OFFSET fft64Kp_2PPRO
	DD			OFFSET fft64Kp_3PPRO, OFFSET fft64Kp_4PPRO
	DD			127*65536+8*256+16, 1
	DD			127*65536+8*256+4, 1
	DD			128*256+16, 1, 1, 1, 0
	DD	2655500, 131072, 0.0719,	1380096
	DD			OFFSET fft128Kp_1, OFFSET fft128Kp_2
	DD			OFFSET fft128Kp_3, OFFSET fft128Kp_4
	DD			OFFSET fft128Kp_1PPRO, OFFSET fft128Kp_2PPRO
	DD			OFFSET fft128Kp_3PPRO, OFFSET fft128Kp_4PPRO
	DD			1, 1
	DD			255*32768, 1
	DD			7*65536+8*256+16, 1, 1, 1, 0
	DD	5250000, 262144, 0.155,		2757376
	DD			OFFSET fft256Kp_1, OFFSET fft256Kp_2
	DD			OFFSET fft256Kp_3, OFFSET fft256Kp_4
	DD			OFFSET fft256Kp_1PPRO, OFFSET fft256Kp_2PPRO
	DD			OFFSET fft256Kp_3PPRO, OFFSET fft256Kp_4PPRO
	DD			1, 1
	DD			511*32768, 1
	DD			15*65536+8*256+16, 2, 1, 1, 0
	DD	10320000, 524288, 0.322,	5514240
	DD			OFFSET fft512Kp_1, OFFSET fft512Kp_2
	DD			OFFSET fft512Kp_3, OFFSET fft512Kp_4
	DD			OFFSET fft512Kp_1PPRO, OFFSET fft512Kp_2PPRO
	DD			OFFSET fft512Kp_3PPRO, OFFSET fft512Kp_4PPRO
	DD			1, 1
	DD			1023*32768, 1
	DD			31*65536+8*256+16, 1, 1, 1, 0
	DD	20400000, 1048576, 0.681,	11023360
	DD			OFFSET fft1024Kp_1, OFFSET fft1024Kp_2
	DD			OFFSET fft1024Kp_3, OFFSET fft1024Kp_4
	DD			OFFSET fft1024Kp_1PPRO, OFFSET fft1024Kp_2PPRO
	DD			OFFSET fft1024Kp_3PPRO, OFFSET fft1024Kp_4PPRO
	DD			1, 1
	DD			2047*32768, 1
	DD			63*65536+8*256+16, 2, 1, 1, 0
	DD	40250000, 2097152, 1.380,	22050816
	DD			OFFSET fft2048Kp_1, OFFSET fft2048Kp_2
	DD			OFFSET fft2048Kp_3, OFFSET fft2048Kp_4
	DD			OFFSET fft2048Kp_1PPRO, OFFSET fft2048Kp_2PPRO
	DD			OFFSET fft2048Kp_3PPRO, OFFSET fft2048Kp_4PPRO
	DD			1, 1
	DD			4095*32768, 1
	DD			127*65536+8*256+16, 1, 1, 1, 0
	DD	79300000, 4194304, 2.919,	44087296
	DD			OFFSET fft4096Kp_1, OFFSET fft4096Kp_2
	DD			OFFSET fft4096Kp_3, OFFSET fft4096Kp_4
	DD			OFFSET fft4096Kp_1PPRO, OFFSET fft4096Kp_2PPRO
	DD			OFFSET fft4096Kp_3PPRO, OFFSET fft4096Kp_4PPRO
	DD			1, 1
	DD			8191*32768, 1
	DD			255*65536+8*256+16, 2, 1, 1, 0
_TEXT32A ENDS

	ASSUME  CS: _TEXT32A, DS: _TEXT32A, SS: _TEXT32A, ES: _TEXT32A

INCLUDE	unravel.mac
INCLUDE	lucas1.mac
INCLUDE pfa.mac
INCLUDE memory.mac
INCLUDE setup.mac
INCLUDE normal.mac

_TEXT32A	SEGMENT

	PUBLIC	_gwinfo1
_gwinfo1 PROC NEAR
	push	esi

; Decide which table to scan

	mov	esi, OFFSET jmptable	; Assume mersenne mod FFTs
	cmp	_INFT, 0		; Check 2^N+1 flag
	jz	short mmod		; Yes, this is 2^N-1 math
	mov	esi, OFFSET jmptablep	; Do 2^N+1 mod FFTs

; Find the table entry using either the specified fft length or
; the exponent being tested.

mmod:	mov	ecx, _INFF		; FFT length to lookup (or zero)
	cmp	ecx, 0			; Was a specific fft length given?
	je	short fexp		; No, search table using exponent
	mov	eax, -1			; Invalidate searching by exponent
	jmp	short flp		; Start searching by fft length
fexp:	mov	eax, _INFP		; Exponent to lookup
	mov	ecx, -1			; Invalidate searching by fft length
flp:	mov	edx, [esi]		; Load maximum exponent
	cmp	eax, edx		; Is our exp less than maximum exp?
	jbe	short fdn		; Yes, we've found our table entry
	mov	edx, [esi+4]		; Load fft length
	cmp	ecx, edx		; Is this the fftlen caller wanted?
	jbe	short fdn		; Yes, we've found our table entry
flp1:	cmp	DWORD PTR [esi], 0	; Look for zero terminator
	lea	esi, [esi+4]		; Point to next word
	jnz	short flp1		; Scan until zero found
	jmp	short flp		; Test next table entry

; We've found our table entry.  Return the address in INFT.

fdn:	mov	_INFT, esi
	pop	esi
	ret
_gwinfo1 ENDP


; Setup routine

	PUBLIC	_gwsetup1
_gwsetup1 PROC NEAR
	push	ebp
	push	ebx
	push	edi
	push	esi

; Initialize

	fninit				; Init FPU

; Compute square root of 0.5

	fld	HALF
	fsqrt
	fstp	SQRTHALF

; Compute the rounding constant

	mov	edi, OFFSET BIGVAL
	mov	DWORD PTR [edi], 12
	fild	DWORD PTR [edi]		;; 3*2^2
	mov	DWORD PTR [edi], 4096
	fimul	DWORD PTR [edi]		;; 3*2^14
	fimul	DWORD PTR [edi]		;; 3*2^26
	fimul	DWORD PTR [edi]		;; 3*2^38
	fimul	DWORD PTR [edi]		;; 3*2^50
	fimul	DWORD PTR [edi]		;; 3*2^62
	fstp	DWORD PTR [edi]

; Compute the constants used in prime factor FFTs

	pfa_5_setup
	pfa_6_setup
	pfa_7_setup

; Find the table entry using either the specified fft length or
; the exponent being tested.

	mov	eax, _PARG		; Exponent we are testing
	mov	_INFP, eax		; Place it where gwinfo looks for it
	mov	eax, _FFTLEN		; Specific fft length to use (or zero)
	mov	_INFF, eax		; Place it where gwinfo looks for it
	mov	eax, _PLUS1		; Type of FFT
	mov	_INFT, eax		; Place it where gwinfo looks for it
	CALL_X	_gwinfo1		; Lookup in table
	mov	esi, _INFT

; We've found our table entry.  Return the amount of memory that needs
; to be allocated.

	mov	eax, [esi+12]

; Return

	pop	esi
	pop	edi
	pop	ebx
	pop	ebp
	ret
_gwsetup1 ENDP


;; Phase 2 of setup.  Caller has allocated memory for us to use.

	PUBLIC	_gwsetup2
_gwsetup2 PROC NEAR
	push	ebp
	push	ebx
	push	edi
	push	esi

; Copy info from the table entry

	mov	esi, _INFT		; Reload table pointer
	mov	ecx, [esi+4]		; Save FFTLEN
	mov	_FFTLEN, ecx
	cmp	_CPU_TYPE, 6		; See if we should run PPRO routines
	je	short ppro		; Yes if a PPRO
	cmp	_CPU_TYPE, 8		; See if we should run PPRO routines
	jge	short ppro		; Yes if a PII, Celeron, PIII, etc
	mov	ecx, [esi+16]		; Save 4 fft routine offsets
	mov	_GWPROCPTRS, ecx	; Forward FFT routine
	mov	ecx, [esi+20]
	mov	_GWPROCPTRS+4, ecx	; Squaring routine
	mov	ecx, [esi+24]
	mov	_GWPROCPTRS+8, ecx	; Multiply (one value already FFTed)
	mov	ecx, [esi+28]
	mov	_GWPROCPTRS+12, ecx	; Multiply (both values already FFTed)
	jmp	short counts		; Now go load the counters
ppro:	mov	ecx, [esi+32]		; Save 4 fft routine offsets
	mov	_GWPROCPTRS, ecx	; Forward FFT routine
	mov	ecx, [esi+36]
	mov	_GWPROCPTRS+4, ecx	; Squaring routine
	mov	ecx, [esi+40]
	mov	_GWPROCPTRS+8, ecx	; Multiply (one value already FFTed)
	mov	ecx, [esi+44]
	mov	_GWPROCPTRS+12, ecx	; Multiply (both values already FFTed)
counts:	mov	ecx, [esi+48]		; Save 2 normalize counters
	mov	addcount1, ecx
	mov	ecx, [esi+52]
	mov	addcount2, ecx
	mov	ecx, [esi+56]		; Save 2 normalize counters
	mov	normcount1, ecx
	mov	ecx, [esi+60]
	mov	normcount2, ecx
	mov	ecx, [esi+64]		; Save 4 counters
	mov	count1, ecx
	mov	ecx, [esi+68]
	mov	count2, ecx
	mov	ecx, [esi+72]
	mov	count3, ecx
	mov	ecx, [esi+76]
	mov	count4, ecx
	lea	esi, [esi+80]		; save table addr for later
	push	esi

; Generate the power-of-two sine/cosine tables used in pass 2 FFTs

	mov	esi, sincos_real	;; Address of sine-cosine array
	sincos_real_setup
	mov	esi, sincos_complex	;; Address of sine-cosine array
	sincos_complex_setup

; Compute two-to-phi and two-to-minus-phi normalization multipliers

	mov	eax, _SRCARG		; Load next available address
	add	eax, 31			; Make it a cache line boundary
	and	eax, 0FFFFFFE0h
	mov	_NORM_ARRAY2, eax	; Store the column data there
	add	eax, 128*NMD		; There are up to 128 columns
	mov	_NORM_ARRAY1, eax	; Store the group data there
	normalize_setup
	mov	_SRCARG, edi		; Save address for next table

; Generate the premultipliers used in pass 1

	mov	esi, _SRCARG		; Load next available address
	add	esi, 31			; Make it a cache line boundary
	and	esi, 0FFFFFFE0h
	mov	pass1_premults, esi	; Store the premult data there
	mov	eax, _FFTLEN		; Load params based on FFT length
	cmp	eax, 65536		; Are there any pass1 multipliers?
	JLE_X	nop1			; No for FFTs below 65536 in size
	mov	ebx, 64			; Assume 64 elements in pass 1
	cmp	count2, 2		; Are there 64 or 128 elements in pass1
	jne	short is64		; Yes, it is 64 elements in pass 1
	mov	ebx, 128		; No, there are 128 elements in pass 1
is64:	shr	eax, 8			; Compute the size of the pass 0/1 FFT
	premultiplier_setup
nop1:	mov	_SRCARG, esi		; Save address for next table

; Generate the premultipliers used in pass 2

	mov	esi, _SRCARG		; Load next available address
	add	esi, 31			; Make it a cache line boundary
	and	esi, 0FFFFFFE0h
	mov	pass2_premults, esi	; Store the premult data there
	mov	eax, _FFTLEN		; Load params based on FFT length
	cmp	eax, 1024		; Are there any pass2 multipliers?
	JLE_X	nop2			; No for FFTs below 1024 in size
	mov	ebx, 64			; Assume 64 elements in pass 2
	cmp	eax, 8192		; Are there 64 or 256 elements in pass2
	jle	short was64		; Yes, it is 64 elements in pass 2
	mov	ebx, 256		; No, there are 256 elements in pass 2
was64:	premultiplier_setup
nop2:	mov	_SRCARG, esi		; Save address for next table

; Compute up to 4 arrays of sin/cos values for PFA FFTs

	mov	edi, OFFSET sincos1	; Address of table pointers
	mov	DWORD PTR [edi], sincos_real
	mov	DWORD PTR [edi+4], sincos_real
	mov	DWORD PTR [edi+8], sincos_real
	mov	DWORD PTR [edi+12], sincos_real
	mov	esi, _SRCARG		; Load next available address
	add	esi, 31			; Make it a cache line boundary
	and	esi, 0FFFFFFE0h
slp:	pop	ebp			; Get pointer to number of PFA elements
	mov	ebx, [ebp]		; Get number of PFA elements
	cmp	ebx, 0
	JZ_X	sdn			; No more data to generate
	lea	ebp, [ebp+4]		; Point to next entry
	mov	[edi], esi		; Save address of this array
	push	ebp
	push	edi
	sc0_setup
	pop	edi
	lea	edi, [edi+4]
	JMP_X	slp			; Look for more to generate
sdn:	mov	_SRCARG, esi		; Save address for next table

; Create multipliers for 2^N+1 arithmetic

	cmp	_PLUS1, 0
	JZ_X	minus1
	mov	esi, _SRCARG		; Load next available address
	add	esi, 31			; Make it a cache line boundary
	and	esi, 0FFFFFFE0h
	mov	plus1_premults, esi
	plus1_mult_setup
	mov	_SRCARG, esi		; Save address for next table
minus1:

; Compute extra bits (the number of adds we can tolerate without
; a normalization operation).  Studies show that exponents below
; maxp - fftlen/2 - 4 can withstand one addition without normalization.

	mov	esi, _INFT
	mov	edx, DWORD PTR [esi]	; Load maximum exponent for fft size
	sub	edx, _PARG
	mov	eax, _FFTLEN
	shr	eax, 1
	add	edx, eax
	sub	edx, 4
	mov	eax, edx
	sub	edx, edx
	div	_FFTLEN
	mov	extra_bits, eax

; Set pointers to add/sub/copy routines

	mov	eax, _FFTLEN
	mov	edi, OFFSET _GWPROCPTRS
	cmp	eax, 128
	jg	short notsm
	mov	DWORD PTR [edi+20], OFFSET gwadd1
	mov	DWORD PTR [edi+24], OFFSET gwaddq1
	mov	DWORD PTR [edi+28], OFFSET gwsub1
	mov	DWORD PTR [edi+32], OFFSET gwsubq1
	mov	DWORD PTR [edi+36], OFFSET gwaddsub1
	mov	DWORD PTR [edi+40], OFFSET gwaddsubq1
	JMP_X	ptrdn
notsm:	cmp	eax, 8192
	jg	short notmed
	mov	DWORD PTR [edi+20], OFFSET gwadd2
	mov	DWORD PTR [edi+24], OFFSET gwaddq2
	mov	DWORD PTR [edi+28], OFFSET gwsub2
	mov	DWORD PTR [edi+32], OFFSET gwsubq2
	mov	DWORD PTR [edi+36], OFFSET gwaddsub2
	mov	DWORD PTR [edi+40], OFFSET gwaddsubq2
	jmp	short ptrdn
notmed:	cmp	eax, 65536
	jg	short notlrg
	mov	DWORD PTR [edi+20], OFFSET gwadd3
	mov	DWORD PTR [edi+24], OFFSET gwaddq3
	mov	DWORD PTR [edi+28], OFFSET gwsub3
	mov	DWORD PTR [edi+32], OFFSET gwsubq3
	mov	DWORD PTR [edi+36], OFFSET gwaddsub3
	mov	DWORD PTR [edi+40], OFFSET gwaddsubq3
	jmp	short ptrdn
notlrg:	mov	DWORD PTR [edi+20], OFFSET gwadd4
	mov	DWORD PTR [edi+24], OFFSET gwaddq4
	mov	DWORD PTR [edi+28], OFFSET gwsub4
	mov	DWORD PTR [edi+32], OFFSET gwsubq4
	mov	DWORD PTR [edi+36], OFFSET gwaddsub4
	mov	DWORD PTR [edi+40], OFFSET gwaddsubq4
ptrdn:	mov	DWORD PTR [edi+16], OFFSET gwcopy

; Return

IFDEF MEM_MEASURE
mov eax, _SRCARG
sub eax, _NORM_ARRAY2
ENDIF

	pop	esi
	pop	edi
	pop	ebx
	pop	ebp
	ret
_gwsetup2 ENDP

;;
;; Copy a number with maximal pipelining
;;

gwcopy	PROC NEAR
	push	esi			; U - Save esi
	mov	esi, _SRCARG		; V - Address of first number
	push	edi			; U - Save edi
	mov	edi, _DESTARG		; V - Address of second number
	mov	eax, [esi-4]		; U - Get needs-normalize counter
	mov	ecx, [esi-8]		; V - Load loop counter (size in bytes)
	mov	[edi-4], eax		; U - Store needs-normalize counter
copylp:	mov	eax, [esi][ecx-32]
	mov	edx, [esi][ecx-28]
	mov	[edi][ecx-32], eax
	mov	[edi][ecx-28], edx
	mov	eax, [esi][ecx-24]
	mov	edx, [esi][ecx-20]
	mov	[edi][ecx-24], eax
	mov	[edi][ecx-20], edx
	mov	eax, [esi][ecx-16]
	mov	edx, [esi][ecx-12]
	mov	[edi][ecx-16], eax
	mov	[edi][ecx-12], edx
	mov	eax, [esi][ecx-8]
	mov	edx, [esi][ecx-4]
	mov	[edi][ecx-8], eax
	mov	[edi][ecx-4], edx
	sub	ecx, 32			; Check loop counter
	jnz	short copylp		; Loop if necessary
	pop	edi
	pop	esi
	ret
gwcopy	ENDP

;
; Utility routine to multiply two numbers and then take a modulo
; (32 bit quantities)
;

	PUBLIC	_emulmod
_emulmod PROC NEAR
	mov	eax, _SRCARG		; Integer1
	mul	_SRC2ARG		; Times integer2
	div	_DESTARG		; Mod integer3
	mov	_DESTARG, edx		; Return the remainder
	ret
_emulmod ENDP

;
; Utility routine to initialize the FPU
;

	PUBLIC	_fpu_init
_fpu_init PROC NEAR
	fninit
	ret
_fpu_init ENDP

;
; Utility routine to read the time stamp counter
;

	PUBLIC	_erdtsc
_erdtsc PROC NEAR
	RDTSC
	mov	_DESTARG, eax		; low 32 bits
	mov	_DEST2ARG, edx		; high 32 bits
	ret
_erdtsc ENDP



IFDEF TIMING2
INCLUDE pfa.mac
INCLUDE lucas1.mac
INCLUDE lucas1p.mac
PUBLIC _timeit
_timeit	PROC NEAR
	push	esi
	push	edi
	push	ebp
	push	ebx
	mov	ebx, 0
	mov	ebp, 0
	mov	esi, _SRCARG		; V - Address of first number
	lea	edi, [esi+1024]
	mov	al, 0
	mov	ecx, 250		; 1000 iterations
	align 16
clp1:	disp four_complex_gpm4_fft_0 8, 16, 32
	lea	esi, [esi+64]
	add	al, 256/4
	jnc	clp1
	lea	esi, [esi-256]
	dec	ecx			; Check loop counter
	jnz	clp1			; Loop if necessary
	pop	ebx
	pop	ebp
	pop	edi
	pop	esi
	ret
_timeit	ENDP
ENDIF


_TEXT32A ENDS
END
