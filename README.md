# Multi-sample-VCF-Case-Control-Filter
A tool to filter variants unique to defined cases in a multi-sample VCF
____________________________________________________________________________________________________________________________________________________________________________________

Use at your own risk. I cannot provide support. All information obtained/inferred with this script is without any implied warranty of fitness for any purpose or use whatsoever.

ABOUT:

This program takes as input a multi-sample VCF (uncompressed or compressed with bgzip) and a two column list of samples and their status (case, control or neutral).  The output is a VCF file of variants found in cases but not in controls.  Additional filtering parameters are:  inclusion of a BED file to filter only selected genomic regions, minimum and maximum depth of coverage limits to consider a variant, the ability to allow multi-allelic variants, the ability to allow missing data in samples, and  control over the percent uniqueness (does a variant have to be 100 percent unique in cases, or can the variant appear less than 100 percent of the time). 

PREREQUISITES:
       
    1. C++ compiler (for example GCC/G++) 
    
INSTALLATION:

This program should work on all systems. Download the .cpp file and compile according to your system. For Linux Ubuntu, compile with g++ (g++ -o casecontrol CaseControlVCFFilter_SingleInput_V2.cpp -std=c++11).

USAGE:

./casecontrol [OPTIONS]

<b>Required arguments:</b>

  --vcf <file>           Input multi-sample VCF file (can be .vcf or .vcf.gz)
  
  --samples <file>       Sample status file (sample_name status per line)
  
  -o, --output <file>    Output VCF file

<b>Optional arguments:</b>

  --bed <file>           BED file with regions (optional, no BED = all regions)
  
  --min-depth <int>      Minimum depth (DP) for case/control samples
  
  --max-depth <int>      Maximum depth (DP) for case/control samples
  
  --multi-allelic        Allow multi-allelic variants (0/2, 1/2, etc.)
  
  -a, --allow-missing    Allow missing data (./. treated as neutral)
  
  --percent-unique <pct> Percent uniqueness threshold (default: 100)
  
  -h, --help             Show this help message

<b>Sample status file format:</b>

Prepare a two column plain text file with the sample name in column 1 and a symbol in column 2 indicating the sample status.  + indicates case, - indicates control and 0 means unknown/ignore.  Samples coded with a zero are not considered when viltering the VCF.  Example file format: 

  sample1 +  
  
  sample2 -  
  
  sample3 0 
