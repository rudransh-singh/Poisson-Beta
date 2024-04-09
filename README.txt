Manual Reference Pages - PoissonBeta V1.0

1. NAME

	PoissonBeta - Program for inferring the kinetics of stochastic gene expression from single-cell RNA-seq data

2. SYNOPSIS

	>>fid_out = PoissonBetaRun('test.txt', 'result.txt', 10000)

3. DESCRIPTION

	PoissonBeta is a program for inferring the kinetics of stochastic gene expression from single-cell RNA-seq data.
	It requires MATLAB with Statistics toolbox.

4. COMMANDS AND OPTIONS
	fid_out = PoissonBetaRun(input_file, output_file, Tmax)

	OPTIONS:
	input_file: input file name for single-cell RNA-seq data
		INPUT FILE FORMAT (delimited by tab, please see 'test.txt' for an example)
		COLUMN1: GENE NAME
		COLUMN2: LENGTH OF TRANSCRIPT
		COLUMN3-N: READ COUNT

	output_file: output file name
	Tmax: a maximum number of Gibbs sampling interations

	OUTPUTS:
	fid_out: file identifier of the output file. If failed, fid_out = -1
		COLUMN1: GENE NAME
		COLUMN2: posterior mean Si (not multiplied by transcript length and size factor)
		COLUMN3: Koni
		COLUMN4: Koffi
		COLUMN5: SKoffi (not multiplied by transcript length and size factor)
		COLUMN6: Exi (not multiplied by transcript length and size factor)
		COLUMN7-M: Pij

5. Precompiled version for Windows
	FILES:
	PoissonBetaRun.exe PoissonBetaRun.prj

	INSTALL:
	Download and install MCRInstaller.exe (version 7.11, 64bit) for MATLAB 2009b Windows.

	Add the MCR directory to the environment variable by opening a command prompt and issuing the DOS command:

	set PATH=<mcr_root>\v711\runtime\win64;%PATH% 

	NOTE: <mcr_root> is the directory where MCR is installed on the target machine.         

	USAGE:
	>PoissonBetaRun test.txt result.txt 10000

	