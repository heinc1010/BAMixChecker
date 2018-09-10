try:
	import sys
	import os
	import re
	import time
	import copy
	import commands
	import argparse
	import numpy as np
	from math import log, e
	from os.path import isfile, join
	from subprocess import Popen,PIPE
	from scipy.stats import entropy
	from multiprocessing import Process, Array, Pool
except ImportError:
	sys.exit("##ERROR: Required python package or packages are not installed.\n \
	 Check and install the package or packages with 'pip install'.\
	 - Required python packages : numpy, scipy.stats, multiprocessing")

class VCFinfo:
	def __init__(self,file_path):
		__slots__ = ( 'dic_idx', 'dic_gt')
		self.dic_idx = {}
		self.dic_gt = {}
		self.file_path = file_path

		fr = open(self.file_path,"r")
		current_chr="chr1"
		count_line=0
		for line in fr:
			count_line += 1
			if line.startswith("#"):
				continue
			lis = line.strip().split("\t")
			if lis[0] != current_chr:
				self.dic_idx[current_chr]=count_line-1
				self.dic_gt[current_chr] = {}
				current_chr=lis[0]

		self.dic_idx[current_chr]=count_line
		self.dic_gt[current_chr] = {}
		fr.close()

def get_file_list(dir_path,user_file_list,file_type="BAM"):
	lis_files_whole = []
	lis_files_ans = []
	if (file_type == "BAM") & (user_file_list != ""):
		fr_file_list = open(user_file_list,"r")
		for line_fl in fr_file_list:
			lis_fl = line_fl.strip().split('\t')
			if len(lis_fl) > 1:
				for fl in lis_fl:
					lis_files_whole.apppend(fl)
				lis_fl_names = [ f.split("/")[-1][:3]+"gvcf" for f in lis_fl ]
				lis_files_ans.apppend(lis_fl_names)
			else:
				lis_files_whole.append(lis_fl[0])

		for file_path in lis_files_whole:
			tmp_file_type = file_path.split(".")[-1].upper()
			if tmp_file_type != file_type:
				print "## ERROR: Input files is needed to be .bam file"
				print "## "+file_path+ " is not .bam"
				exit()
	else:
		lis_files_in_dir = [dir_path+f for f in os.listdir(dir_path) if isfile("/".join([dir_path, f]))]
		for file_path in lis_files_in_dir:
			tmp_file_type = file_path.strip().split(".")[-1].upper()
			if tmp_file_type == file_type:
				lis_files_whole.append(file_path)
	if len(lis_files_whole)%2 != 0:
		print "## ERROR: Total number of files is not even - "+str(len(lis_files_whole))
		for f in lis_files_whole:
			print f
		exit()


	return file_type,lis_files_whole,lis_files_ans

def make_bed( total_SNP_bed, targeted_bed , OUTDIR):
	bedtools_path = ""
	fr_config = open(os.path.dirname(os.path.realpath(__file__))+"/MuBaMer.config","r")
	for line in fr_config:
		if line.startswith("BEDTOOLS="):
			bedtools_path = line.strip().split("=")[1]
			break
	fr_config.close()
	if bedtools_path == "":
		print "##ERROR: BEDTOOLS path is not set in 'MuBaMer.config' file."
		exit()
	os.system("{0} intersect -b {1} -a {2} > {3}targetSNPs.bed".format(bedtools_path,targeted_bed,total_SNP_bed,OUTDIR))
	(exitstatus, line_count) = commands.getstatusoutput("cat {0}targetSNPs.bed | wc -l".format(OUTDIR))
	if int(line_count) < 200:
		if total_SNP_bed not in ["gnomad_hg38_AF10.bed","gnomad_hg19_AF10.bed"]:
			os.system("rm {0}targetSNPs.bed".format(OUTDIR))
			return None
		else:
			print "WARNING: The number of target SNP sites to compare is under 200."
			print "*Make a custermized list of target SNPs - targetSNPs.bed\n",
			return OUTDIR + "targetSNPs.bed"
	else:
		print "*Make a custermized list of target SNPs - targetSNPs.bed\n",
#		print total_SNP_bed
		return OUTDIR + "targetSNPs.bed"

def get_gt(f1):
	fr = open(f1.file_path,"r")
	for line in fr:
		if line.startswith("#"):
			continue
		dp_idx = -1
		pl_idx = -1
		lis = line.strip().split("\t")
		lis_fm = lis[8].strip().split(":")
		for i in range(0,len(lis_fm)):
			if lis_fm[i] == "DP":
				dp_idx = i
			elif lis_fm[i] == "PL":
				pl_idx = i
		lis_info = lis[9].strip().split(":")
		if int(lis_info[dp_idx]) < 5:
			continue
		lis_pl = [ int(pl) for pl in lis_info[pl_idx].split(",") ]
		if len(lis_pl) > 6:
			continue

		lis_pl = lis_pl[:3]
		baf = int(lis_pl.index(min(lis_pl))*5.0)
		f1.dic_gt[lis[0]][lis[1]] = baf
	fr.close()

def cal_cor_each(f1,f2,cor_arr,x):
	count_all=0
	count_match=0
	count_diff_baf = 0
	lis_diff_snp = []
	for ch in f1.dic_gt.keys():
		for pos in f1.dic_gt[ch].keys():
			baf1 = f1.dic_gt[ch][pos]
			try:
				baf2 = f2.dic_gt[ch][pos]
			except:
				continue
			count_all += 1
			if baf1 == baf2:
				count_match += 1
			else:
				count_diff_baf += 1
				lis_diff_snp.append(ch+":"+pos)

	if count_match == 0:
		fin_cor = 0
	else:
		fin_cor = round(count_match*1.0/count_all,4)
	cor_arr[x] = fin_cor

def cal_cor(lis_files, max_prc):
	cor_matrix = []
	lis_vcf = []
	for i in range(0,len(lis_files)):
		lis_vcf.append(VCFinfo(lis_files[i]))

	procs = []
	for i in range(0,len(lis_vcf)):
		get_gt(lis_vcf[i])

	for i in range(0,len(lis_vcf)):
		cor_matrix.append([])
		procs = []
		cor_arr = Array('d',[-1]*len(lis_vcf))
		for j in range(0, len(lis_vcf)):
			if j < i :
				cor_matrix[i].append(cor_matrix[j][i])
				continue
			elif j > i :
				procs.append(Process(target=cal_cor_each, args=(lis_vcf[i],lis_vcf[j],cor_arr,j)))
				if len(procs) == max_prc:
					for p in procs:
						p.start()
					for p in procs:
						p.join()
					procs=[]
			else:
				cor_matrix[i].append(-1)
				continue
		if procs != []:
			for p in procs:
				p.start()
			for p in procs:
				p.join()
		for cor in cor_arr[i+1:]:
			cor_matrix[i].append(cor)

	return cor_matrix

def pairing(cor_matrix, lis_files):
	num_files = len(lis_files)
	lis_m_idx =[]
	set_skip = set([])
	smp_pairs = {}
	for i in range(0,len(lis_files)):
		smp_pairs[lis_files[i]] = [ lis_files[j] for j in range(0,len(lis_files)) if cor_matrix[i][j] > 0.7 ]
	set_tmp_g = []
	for f1 in smp_pairs.keys():
		lis_tmp = [f1]
		for f2 in smp_pairs[f1]:
			lis_tmp.append(f2)
		lis_tmp.sort()
		set_tmp_g.append(lis_tmp)
	lst_smp_pairs = {}
	for g in set_tmp_g:
		lst_smp_pairs[g[0]] = g[1:]

	return lst_smp_pairs

def compareVectors(v1, v2, entropies):
	score = 0
	for i in range(0, len(v1)):
		if v1[i]==v2[i]:
			score += entropies[i]
		else:
			score -= entropies[i]
	return score

def findmatch(scorevector):
	maxitem = scorevector[0][0]
	maxscore = scorevector[0][1]
	for i in range(1, len(scorevector)):
		if scorevector[i][1] > maxscore:
			maxitem = i
			maxscore = scorevector[i][1]
	return maxitem, maxscore

def entropy1(labels, base=None):
	value, counts = np.unique(labels, return_counts=True)
	return entropy(counts, base=base)

def get_max(lis_scores):
	lis_vs = lis_scores[:]
	maxscore = lis_vs[0][1]
	lis_t_topscore = []
	for t_vs in lis_vs:
		if t_vs[1]>=maxscore:
			maxscore = t_vs[1]
	for t_vs in lis_vs:
		if t_vs[1]==maxscore:
			lis_t_topscore.append(t_vs[0])
	return lis_t_topscore

def get_sw_pairs(lis_vcf_files,smp_pairs):
	separators = r'_|-|\.'
	dic_split_samples = {}
	dic_sample_scores = {}
	dic_topscores = {}
	lis_factorsize = []
	for v in lis_vcf_files:
		factors = re.split(separators, v)
		lis_factorsize.append(len(factors))
		dic_split_samples[v] = factors
	if len(set(lis_factorsize)) != 1:
		return None,None
	len_fct = lis_factorsize[0]
	lis_factors = []
	lis_entropies = []
	for i in range(0, len_fct):
		lis_factors.append([])
	for v in dic_split_samples.keys():
		fv = dic_split_samples[v]
		for i in range(0, len_fct):
			lis_factors[i].append(fv[i])
	for i in range(0, len_fct):
		lis_entropies.append(entropy1(lis_factors[i]))
	for s1 in dic_split_samples.keys():
		dic_sample_scores[s1] = []
		for s2 in dic_split_samples.keys():
			if s1!=s2:
				score = compareVectors(dic_split_samples[s1], dic_split_samples[s2], lis_entropies)
				dic_sample_scores[s1].append((s2, score))
	dic_sw = {}
	dic_un_p = {}
	for f1 in smp_pairs.keys(): 
		if smp_pairs[f1] == []:
			dic_un_p[f1] = get_max(dic_sample_scores[f1])
		else:
			dic_sw[f1] = []
			lis_m_f = get_max(dic_sample_scores[f1])
			for f2 in smp_pairs[f1]:
				if f2 not in lis_m_f:
					dic_sw[f1].append(f2)
	return dic_sw, dic_un_p

def get_sw_pairs_ans(lis_vcf_files, smp_pairs, lis_ans):
	f1s = smp_pairs.keys()
	f1s.sort()
	if f1 in f1s:
		for g in lis_ans:
			if f1 in g:
				g_f = copy.deepcopy(g)
				g_f.remove(f1)
				g_f.sort()
				smp_pairs[f1].sort()
				if g_f != smp_pairs[f1]:
					if smp_pairs[f1] == []:
						dic_un_p[f1] = g_f
					else:
						for f2 in g_f:
							dic_sw[f1].append(f2)
	return dic_sw, dic_un_p

def make_result_file_no_file_name_info(cor_matrix,smp_pairs,lis_vcf_files,OutputDIR,lis_ans,lis_files):
	print "          Skip making 'Unmatched_pair.txt' file"
	fw_m_m = open(OutputDIR+"Matched_pair.txt","w")
	fw_m_m.write("#Matched pair by genotype.\n")
	for i in range(0,len(lis_files)):
		f1 = smp_pairs.keys()[i]
		if smp_pairs[f1]  == []:
			continue
		for f2 in smp_pairs[f1]:
			if i <= lis_vcf_files.index(f2):
				continue
			flag_less_informative = False
			fw_m_m.write(f1+"\t"+f2+"\t")
			score = cor_matrix[lis_vcf_files.index(f1)][lis_vcf_files.index(f2)]
			if score < 0.8:
				flag_less_informative = True
			fw_m_m.write(str(score)+"\n")
			if flag_less_informative:
				fw_m_m.write("-> *This pair scores under 0.8 which is less informative. The 'less informative score' dosen't mean that the pair is not matched but may have some problem of purity or copy number etc. in the sample.\n")
	fw_m_m.close()

def make_result_file(cor_matrix,smp_pairs,lis_vcf_files,OutputDIR,lis_ans):
	return_v = 1
	fw_a_m = open(OutputDIR+"Total_result_matrix.txt","w")
	lis_files = smp_pairs.keys()
	lis_files.sort()
	if ( lis_ans == [] ) & (len(lis_vcf_files) <= 6) :
		print "Warning : The number of files is not enough to pair by file names."
		make_result_file_no_file_name_info(cor_matrix,smp_pairs,lis_vcf_files,OutputDIR,lis_ans,lis_files)
		return_v = 0
	else:
		count_m = 0
		count_s = 0
		count_u = 0
		if lis_ans == []:
			dic_sw, dic_un_p = get_sw_pairs(lis_vcf_files,smp_pairs)
			if ( dic_sw == None ) & ( dic_un_p == None):
				print "Warning : The file names don't have detectable regulation."
				make_result_file_no_file_name_info(cor_matrix,smp_pairs,lis_vcf_files,OutputDIR,lis_ans,lis_files)
				return_v = 0
				len_v = len(lis_vcf_files)
				for i in range(0,len_v-1):
					for j in range(i+1, len_v):
						m_um="unmatched"
						if cor_matrix[i][j] > 0.7:
							m_um = "matched"
						fw_a_m.write(lis_vcf_files[i]+"\t"+lis_vcf_files[j]+"\t"+str(cor_matrix[i][j])+"\t"+m_um+"\n")
				fw_a_m.close()
				return return_v
			else:
				print "Dtected pairs by file names."
		else:
			dic_sw, dic_un_p = get_sw_pairs_ans(lis_vcf_files, smp_pairs, lis_ans)
		fw_m_m = open(OutputDIR+"Matched_pair.txt","w")
		fw_s_m = open(OutputDIR+"Unmatched_pair.txt","w")
		perfect_m = True
		for f1 in lis_files:
			if len(smp_pairs[f1]) > 0:
				for sw in dic_sw[f1]:
					perfect_m = False
					break
			else:
				for un in dic_un_p:
					perfect_m = False
					break
		if perfect_m:
			fw_s_m.write("No swapped or unpaired file.")
		else:
			return_v = 2
			fw_s_m.write("#List of pairs are not matched by name but by genotype.\n")
		fw_m_m.write("#Matched pair by both genotype and name.\n")
		for i in range(0,len(lis_files)):
			f1 = lis_files[i]
			if smp_pairs[f1]  == []:
				continue
			if len(smp_pairs[f1]) > 1:
				print smp_pairs[f1]
			for f2 in smp_pairs[f1]:
				flag_less_informative = False
				if f2 in dic_sw[f1]:
					count_s +=  1
					fw_s_m.write(f1+"\t"+f2+"\t")
					score = cor_matrix[lis_vcf_files.index(f1)][lis_vcf_files.index(f2)]
					if score < 0.8:
						flag_less_informative = True
					fw_s_m.write(str(score)+"\n")
					if flag_less_informative:
						fw_s_m.write("-> *This pair scores under 0.8 which is less informative. The 'less informative score' dosen't mean that the pair is not matched but may have some problem of purity or copy number etc. in the sample.\n")
				else:
					count_m += 1
					flag_less_informative = False
					fw_m_m.write(f1+"\t"+f2+"\t")
					score = cor_matrix[lis_vcf_files.index(f1)][lis_vcf_files.index(f2)]
					if score < 0.8:
						flag_less_informative = True
					fw_m_m.write(str(score)+"\n")
					if flag_less_informative:
						fw_m_m.write("-> *This pair scores under 0.8 which is less informative. The 'less informative score' dosen't mean that the pair is not matched but may have some problem of purity or copy number etc. in the sample.\n")
		if not perfect_m:
			fw_s_m.write("\n#List of samples are matched with nothing by genotype.\n")
			for i in range(0,len(lis_files)):
				f1 = lis_files[i]
				if smp_pairs[f1]  == []:
					count_u += 1
					fw_s_m.write(f1+"\n")
					for f2 in dic_un_p[f1]:
						fw_s_m.write("-> pair by name with "+ f2 +" ( score : "+str(cor_matrix[lis_vcf_files.index(f1)][lis_vcf_files.index(f2)]) +")\n")
		fw_s_m.close()
		fw_a_m.write("#Matched pairs : {0}, Swaped pairs : {1}, Unpaired sample : {2}\n".format(count_m, count_s,count_u))
		fw_m_m.close()
	len_v = len(lis_vcf_files)
	for i in range(0,len_v-1):
		for j in range(i+1, len_v):
			m_um="unmatched"
			if cor_matrix[i][j] > 0.7:
				m_um = "matched"
			fw_a_m.write(lis_vcf_files[i]+"\t"+lis_vcf_files[j]+"\t"+str(cor_matrix[i][j])+"\t"+m_um+"\n")
	fw_a_m.close()

	return return_v

def run_p(cmm):
	print cmm
	try:
		prc= Popen(cmm, stdout=PIPE, shell=True, stderr=PIPE)
		stdoutput, stderr = prc.communicate()
		if prc.returncode != 0:
			print stderr
			exit()
	except KeyboardInterrupt:
		for p in multiprocessing.active_children():
			p.terminate()
		exit()

def run_bam(lis_bam_files,OutputDIR,reference_file,bed_file,max_prc,lis_ans):
	print "*Max process : ",
	print max_prc
	max_prc = int(max_prc)
	HC_path = ""
	fr_config = open(os.path.dirname(os.path.realpath(__file__))+"/MuBaMer.config","r")
	for line in fr_config:
		if line.startswith("GATK="):
			HC_path = line.strip().split("=")[1]
			break
	fr_config.close()

	if HC_path == "":
		print "##ERROR: GATK path is not set in 'MuBaMer.config' file."
		exit()

	vcf_file_path = OutputDIR + "HaplotypeCaller/"
	cmm = "mkdir -p {0}".format(vcf_file_path)
	os.system(cmm)

	lis_cmm = []
	for i in range(0,len(lis_bam_files)):
		vcf_file = "{0}{1}gvcf".format(vcf_file_path, lis_bam_files[i].split("/")[-1][:-3])
		cmm ="{0} HaplotypeCaller -I {1} -O {2} -R {3} -L {4} -ERC BP_RESOLUTION".format(HC_path,lis_bam_files[i],vcf_file,reference_file,bed_file)
		lis_cmm.append(cmm)
	pool = Pool(processes = max_prc)
	try:
		pool.map(run_p,lis_cmm)
		pool.close()
		pool.join()
	except KeyboardInterrupt:
		print "##ERROR: KeyboardInterrupt"
		pool.close()
		pool.terminate()
		pool.join()
		exit()

	lis_vcf_files = []
	for f in lis_bam_files:
		tmp_f = f.split("/")[-1]
		tmp_f = tmp_f[:-3]+"gvcf"
		lis_vcf_files.append(vcf_file_path+tmp_f)

	lis_vcf_files.sort()
	cor_matrix = cal_cor(lis_vcf_files, max_prc)

	for i in range(0,len(lis_vcf_files)):
		lis_vcf_files[i] = lis_vcf_files[i].split("/")[-1]

	smp_pairs = pairing(cor_matrix,lis_vcf_files)
	result = make_result_file(cor_matrix,smp_pairs,lis_vcf_files,OutputDIR,lis_ans)
	if result == 1:
		print "Perfect match."
	elif result == 2:
		print "Swapped file exist. Check a 'swapped_pair.txt'."


if __name__ == "__main__":
	start_t = time.time()

	parser = argparse.ArgumentParser(prog="MuBaMer", description="Bam file swap checker. Grouping files from same indivisual.")
	parser.add_argument('-d','--DIR', default="", help="Directory path of the .BAM files")
	parser.add_argument('-l', '--List', default="", help="A file with the list of files")
	parser.add_argument('-r', '--Ref', default="",required=True, help="Reference file")
	parser.add_argument('-o','--OutputDIR', default="", help="Output directory path")
	parser.add_argument('-b','--BEDfile', default="", help=".bed file for Targeted sequencing data.")
	parser.add_argument('-v', '--RefVer', default="hg38", choices=['hg38','hg19'], help="Version of reference : 'hg19' or 'hg38'. Default = 'hg38'")
	parser.add_argument('-p', '--MaxProcess', default="5", help="The number of max process. Default = 5")
	parser.add_argument('--RemoveVCF',action='store_true', help="Use to remove called germline VCF during sample identification in bam file mode after running.")

	args = parser.parse_args()
	if args == None:
		print "##ERROR: -h ot --help option to get an information to use MuBaMer."
		exit()
	if args.DIR == "":
		if args.List == "":
				print "##ERROR: There is no information about input files. Use -d or -l option for input file information."
				exit() ## error code
	else:
		if args.List != "":
			print "##ERROR: Option -d and -l are excusive. Try with one of the options."
			exit() ## error code
	if not args.DIR.endswith("/"):
		args.DIR = args.DIR+"/"

	if args.OutputDIR == "":
		args.OutputDIR = args.DIR+"/MuBaMer/"
	elif not args.OutputDIR.endswith("/"):
		args.OutputDIR += "/"

	cmm = "mkdir -p {0}".format(args.OutputDIR)
	os.system(cmm)

	if args.Ref == "":
		print "ERROR: Reference file is necessary for BAM file mode. Use -r option."
		print exit()

	flag_chr = False
	(exitstatus, header) = commands.getstatusoutput("head -1 {0}".format(args.Ref))
	if header.startswith(">chr"):
		flag_chr = True

	if flag_chr:
		bed_file_path = os.path.dirname(os.path.realpath(__file__))+"/bed/"
	else:
		bed_file_path = os.path.dirname(os.path.realpath(__file__))+"/bed/noChr/"

	bed_file = None
	if args.RefVer in ["hg38","hg19"]:
		if args.BEDfile != "":
			print "*Run for targeted sequecing data"
			bed_file = make_bed("{0}gnomad_{1}_AF{2}_AF{3}_All.bed".format(bed_file_path,args.RefVer,45,35), args.BEDfile , args.OutputDIR)
			for AF in range(45,5,-5):
				for AF_all in range(AF,0,-10):
					if AF_all >= AF:
						continue
					AF_all = int(AF_all/10)*10
					if AF_all != 0:
						bed_file = make_bed("{0}gnomad_{1}_AF{2}_AF{3}_All.bed".format(bed_file_path,args.RefVer,AF,AF_all), args.BEDfile , args.OutputDIR)
					else:
						bed_file = make_bed("{0}gnomad_{1}_AF{2}.bed".format(bed_file_path,args.RefVer,AF), args.BEDfile , args.OutputDIR)
					if bed_file != None:
						break
				if bed_file != None:
					break
		else:
			bed_file = "{0}gnomad_{1}_AF45_AF35_All.bed".foamt(bed_file_path,args.RefVer)
	else:
		print "##ERROR: Option -v should be 'hg19' or 'hg38'."
		exit()

	file_type,lis_files,lis_ans = get_file_list(args.DIR, args.List)
	run_bam(lis_files,args.OutputDIR,args.Ref,bed_file,args.MaxProcess,lis_ans)

	if args.RemoveVCF:
		cmm = "rm -r {0}HaplotypeCaller/".format(args.OutputDIR)
		os.system(cmm)

	end_t = time.time() -  start_t
	print "Running time: "+str(round(end_t/60,2))+" min"

