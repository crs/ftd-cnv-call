# read all datasets filenames from folder

STUDY = 'study.list.660only'

with open(STUDY) as f:
	FILES = f.readlines()
FILES = [x.strip() for x in FILES]

#FILES, = glob_wildcards("PennCNV-Exports/{dataset}.txt.gz")
#DATASETS, SAMPLES = glob_wildcards("processing/{dataset}/Sample.{sample}")
PENNCNV = '/home/chris/CNV/penncnv/'
PENNCNVLIB = PENNCNV + 'lib/'
GCMODEL, HMM, PFB = [PENNCNVLIB + file for file in ['hhall.hg18.gcmodel', 'hhall.hmm', 'hhall.hg18.pfb']]


FILTERCNV = 'perl ' + PENNCNV + 'filter_cnv.pl '
SCAN = 'perl ' + PENNCNV + 'scan_region.pl '
TOPLINK = 'perl ' + PENNCNV + 'penncnv_to_plink.pl '
VISUALIZECNV = 'perl ' + PENNCNV + 'visualize_cnv.pl '

#GAPFILE= '/home/chris/CNV/dbGaP/centromers.telomers.txt'
GAPFILE= '/home/chris/CNV/dbGaP/allGaps.txt' 

REF = '/home/chris/CNV/ref/'

REFGENE = REF + 'refGene.txt.hg18'
REFLINK = REF + 'refLink.txt.hg18'
HG= REF + 'glist-hg18_mod'

CASEFAM='FTD/FTD_660_sort.fam'
CASECNV='FTD/FTD_660_sort.cnv'

MINSNP = 10
MINLEN = '100k'
CONFIDENCE = 20.0



#stringent thresholds
RRSD=0.235
BAFD=0.01
WF=0.05
QCNUMCNV = 100

#RRSD=0.3
#BAFD=0.002
#WF=0.04


def getSamplesForDataset(wildcards):
	filename = 'unpacked/'+wildcards.dataset + '.txt'
	with open(filename) as f: line = f.readline()
	fields = line.split()
	reg = re.compile(r"^[^.]*")
	return [re.search(reg, f).group(0) for f in fields if 'GType' in f]

def countColumns(filename):
	with file(filename) as f: line = f.readline()
	return len(line.split())

import glob

def getSampleNames(wildcards):
	s, = glob_wildcards('processing/'+ wildcards.dataset + '.{sample}')
	return s
	
def getSampleNames2(wildcards):
	s, = glob_wildcards('processing/'+ wildcards.dataset + '.{sample}')
	return expand('filtered_cnvs/' + wildcards.dataset + '.{sample}.cnv', sample=s)
	
def logfiles(wildcards):
	s, = glob_wildcards('processing/'+ wildcards.dataset + '.{sample}')
	return expand('rawcnv/' + wildcards.dataset + '.{sample}.log', sample=s)

def datasets():
	#samples = glob_wildcards('filtered_cnvs/' + wildcards.dataset + '.{sample}.cnv')
	#print(samples)
	#return glob.glob('filtered_cnvs/'+ wildcards + '.*.cnv')
	print ("DAFUQ")
	return glob_wildcards("processing/{dataset}.*")
	
#def samples(wildcards):
#	samples, = glob_wildcards('filtered_cnvs/' + wildcards.dataset + ".{sample}.cnv")
#	return samples
	
DATASETS,SAMPLES = glob_wildcards("processing/{dataset}.{samples}")
#SAMPLES = glob_wildcards("processing/*.{samples}")
	
#def logfiles(wildcards):
#	return glob.glob('rawcnv/' + wildcards.dataset + '.*.log')

#wildcard_constraints:
#dataset = "|".join(FILES)#,sample = "\w\-\_*"


rule all:	
	input: 	expand('rawcnv/{dataset}.{sample}.rawcnv', zip, dataset=DATASETS,sample=SAMPLES), #dynamic(expand('processing/{dataset}.{{sample}}',dataset=FILES)), #expand('unpacked/{dataset}.txt', dataset=FILES)
			expand('filtered_cnvs/{dataset}.{sample}.cnv', zip, dataset=DATASETS, sample=SAMPLES),
			expand('all_rawcnv/{dataset}.cnv', dataset=FILES),
			expand('all_rawcnv/{dataset}.log', dataset=FILES),
			#expand('qc/{dataset}.rg18', dataset=FILES),
			expand('result/{dataset}.{ext}',dataset=FILES,ext=['fam','cnv.map']),
			expand('final/FTD.{ext}', dataset=FILES, ext=['fam','cnv','rg18','bed']),
			'final/FTD.rawcnv'#,
			#expand('qc/{dataset}.gapsremoved', dataset=DATASETS)
			#expand(dynamic('rawcnv/{dataset}.{{sample}}.rawcnv'),dataset=FILES)

rule createUCSCTrack:
	input: 'final/FTD.rg18'
	output: 'final/FTD.bed'
	shell: VISUALIZECNV +"{input} --format bed --track 'CNV FTD cases' > {output}"
			
rule concatAllStudies:
	input: expand('result/{dataset}.{{ext}}', dataset=FILES)
	output: 'final/FTD.{ext}'
	shell: "cat {input} | sed -e 's,processing/PennCNVexport_all_,,g' > {output}"

rule concatRG18:
	input: expand('qc/{dataset}.rg18', dataset=FILES)
	output: 'final/FTD.rg18'
	shell: "cat {input} | sed -e 's,processing/PennCNVexport_all_,,g' > {output}"
	
rule concatrawcnvs:
	input: expand('qc/{dataset}.rg18', dataset=FILES)
	output: 'final/FTD.rawcnv'
	shell: "cat {input} > {output}"
	
rule intersect:
  input: cnv='result/{dataset}.cnv', map='result/{dataset}.cnv.map', fam='result/{dataset}.fam'
  output: 'result/{dataset}.cnv.indiv', 'result/{dataset}.cnv.regional.summary', 'result/{dataset}.cnv.regional.summary.mperm'
  params: infile='result/{dataset}', outfile='result/{dataset}'
  shell: 'plink --cfile {params.infile} --allow-no-sex --cnv-intersect ' + HG + ' --cnv-test-region --mperm 1 --out {params.outfile}' 
			
rule createFAM:
  input: 'result/{dataset}.cnv'
  output: 'result/{dataset}.fam'
  message: "Create FAM"
  shell: './createCaseFAM.sh {input} {output}'
  
rule makeMap:
  input: 'result/{dataset}.cnv'
  output: 'result/{dataset}.cnv.map'
  params: prefix='result/{dataset}'  
  message: "Create mapfile"
  shell: 'plink --cnv-list {input} --cnv-make-map --out {params.prefix} --allow-no-sex'
			
			
rule convert2plink:
  input: 'qc/{dataset}.gapsremoved'
  output: 'result/{dataset}.cnv'
  message: "Convert PennCNV file {input} to plink compatible format\nOutput in {output}\n"
  shell: TOPLINK + ' -i {input} -o {output}'

 
rule scan:
  input: cnv="qc/{dataset}.gapsremoved"
  output: "qc/{dataset}.rg18"
  shell: SCAN + '{input} -refgene '+ REFGENE + ' -reflink ' + REFLINK + ' --expandmax 5m > {output}'

rule excludeGaps:
	input: gaps='qc/{dataset}.imm', cnv='qc/{dataset}.goodcnv'
	output: "qc/{dataset}.gapsremoved"
	shell: "fgrep -v -f {input.gaps} {input.cnv} > {output}"
  
rule findGaps:
	input: cnv="qc/{dataset}.goodcnv", gap=GAPFILE
	output: imm='qc/{dataset}.imm'
	#shell: SCAN + '{input.cnv} {GAPFILE} -minqueryfrac 0.5 > {output.imm}; fgrep -v -f {output.imm} {input} > {output.goodcnv}'
	shell: SCAN + '{input.cnv} {GAPFILE} --minqueryfrac 0.5 > {output.imm}'
  
rule qc:
  input: cnv='all_rawcnv/{dataset}.cnv', log='all_rawcnv/{dataset}.log'
  output: cnvs='qc/{dataset}.goodcnv', qcpass='qc/{dataset}.qcpass', qcsumout='qc/{dataset}.qcsum'
  params: bafd=BAFD, wf=WF, rrsd=RRSD, qcnumcnv=QCNUMCNV
  shell: FILTERCNV + ' {input.cnv} -qclogfile {input.log} -qclrrsd {params.rrsd} -qcbafdrift {params.bafd} -qcwf {params.wf} -qcnumcnv {params.qcnumcnv} -qcpassout {output.qcpass} -qcsumout {output.qcsumout} -out {output.cnvs}' 
			
			
rule mergeLogs:
  input: logfiles
  output: logfile='all_rawcnv/{dataset}.log'
  shell: 'cat {input} > {output}'
 
 
rule mergeSamples:
  #input: expand('filtered_cnvs/{{dataset}}.{sample}.cnv', sample=getSampleNames)
  input: getSampleNames2
  output: 'all_rawcnv/{dataset}.cnv'
  shell: 'cat {input} > {output}'
  
rule filterCNVs:
  input: 'rawcnv/{dataset}.{sample}.rawcnv'
  output: 'filtered_cnvs/{dataset}.{sample}.cnv'
  params: minsnp=MINSNP, minlen=MINLEN, conf=CONFIDENCE, numcnv=QCNUMCNV
  shell: FILTERCNV + ' --numsnp {params.minsnp} --length {params.minlen} --confidence {params.conf} {input} > {output}'
			
rule call_cnvs:
	input: 'processing/{dataset}.{sample}'
	output: rawcnv='rawcnv/{dataset}.{sample}.rawcnv', logfile='rawcnv/{dataset}.{sample}.log'
	params: detect_cnv='perl '+ PENNCNV + 'detect_cnv.pl', hmm=HMM, gcmodel=GCMODEL, pfb=PFB
	shell: "{params.detect_cnv} --test --hmmfile {params.hmm} --pfbfile {params.pfb} --gcmodelfile {params.gcmodel} --conf --log {output.logfile} --out {output.rawcnv} {input} " #>/dev/null 2>/dev/null

#rule split_samples:
#	input: dataset='unpacked/{dataset}.txt'
#	output: dynamic('processing/{dataset}.{sample}')
#	#output: expand('processing/{{dataset}}.{sample}', sample=getSamplesForDataset)
#	params: kcolumn=PENNCNV+ 'kcolumn.pl', beforestring=".GType", folder='processing'
#	shell: "{params.kcolumn} {input} split 3 -heading 3 -tab -name -beforestring {params.beforestring} -out {params.folder}/{wildcards.dataset}" #>/dev/null 2>/dev/null"""
	
	
#rule unpack_raw_data:
#	input: 'PennCNV-Exports/{dataset}.txt.gz'
#	output: 'unpacked/{dataset}.txt'
#	shell: 'gunzip {input} -c > {output}'