#!/usr/bin/env python2.7

#BAMixChecker version 1.0
#Author : Hein Chun (heinc1010@gmail.com)

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
	from os.path import isfile, join, abspath
	from subprocess import Popen,PIPE
	from scipy.stats import entropy
	from multiprocessing import Process, Array, Pool
except ImportError:
	sys.exit("## ERROR: Required python package or packages are not installed.\n \
	 Check and install the package or packages with 'pip install'.\
	 - Required python packages : numpy, scipy.stats, multiprocessing")

class VCFinfo:
	def __init__(self,file_path):
		__slots__ = ( 'dic_idx', 'dic_gt', 'file_path')
		self.dic_idx = {}
		self.dic_gt = {}
		self.file_path = file_path

		fr = open(self.file_path,'r')
		current_chr='chr1'
		count_line=0
		for line in fr:
			count_line += 1
			if line.startswith('#'):
				continue
			lis = line.strip().split('\t')
			if lis[0] != current_chr:
				self.dic_idx[current_chr]=count_line-1
				self.dic_gt[current_chr] = {}
				current_chr=lis[0]

		self.dic_idx[current_chr]=count_line
		self.dic_gt[current_chr] = {}
		fr.close()

def get_file_list(dir_path,user_file_list):
	lis_files_whole = []
	lis_files_ans = []
	if user_file_list != '':
		fr_file_list = open(user_file_list,'r')
		for line_fl in fr_file_list:
			lis_fl = line_fl.strip().split('\t')
			if len(lis_fl) == 2:
				for fl in lis_fl:
					lis_files_whole.append(fl)
				lis_fl_pair = [ f for f in lis_fl ]
				lis_files_ans.append(lis_fl_pair)
			else:
				lis_files_whole.append(lis_fl[0])

		for file_path in lis_files_whole:
			tmp_file_type = file_path.split(".")[-1].upper()
			if tmp_file_type != 'BAM':
				print "## ERROR: Input files in the list are needed to be .bam file"
				print "## "+file_path+ " is not .bam"
				exit()
	else:
		lis_files_in_dir = [ dir_path+f for f in os.listdir(dir_path) if isfile("/".join([dir_path, f])) ]
		for file_path in lis_files_in_dir:
			tmp_file_type = file_path.strip().split(".")[-1].upper()
			if tmp_file_type == 'BAM':
				lis_files_whole.append(file_path)
	if len(lis_files_whole)%2 != 0:
		print "## ERROR: Total number of files is not even - "+str(len(lis_files_whole))
		for f in lis_files_whole:
			print f
		exit()
	lis_files_whole.sort()
	lis_files_ans.sort()
	return lis_files_whole,lis_files_ans

def make_bed( total_SNP_bed, targeted_bed , OUTDIR, bedtools_path ):
	os.system("{0} intersect -b {1} -a {2} | uniq > {3}targetSNPs.bed".format(bedtools_path,targeted_bed,total_SNP_bed,OUTDIR))
	(exitstatus, line_count) = commands.getstatusoutput("cat {0}targetSNPs.bed | wc -l".format(OUTDIR))
	if int(line_count) < 200:
		if total_SNP_bed.split("/")[-1] not in ["gnomad_hg38_AF10.bed","gnomad_hg19_AF10.bed"]:
			os.system("rm {0}targetSNPs.bed".format(OUTDIR))
			return None
		else:
			print "## WARNING: The target size is too small, so the number of SNP sites to compare is under 200."
			print "##          The specificity could be lower with the small numbr of SNP loci."
			print "Make an optimized list of SNPs to compare - 'targetSNPs.bed'\n",
			return OUTDIR + "targetSNPs.bed"
	else:
		print "Make an optimized list of SNPs to compare - 'targetSNPs.bed'\n",
		return OUTDIR + "targetSNPs.bed"

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

def run_HC( lis_bam_files, OutputDIR, reference_file, bed_file, max_prc, HC_path ):
	print "Calling the variants information with GATK HaplotypeCaller"
	print "Max process : ",
	print max_prc
	max_prc = int(max_prc)
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
		print "## ERROR: KeyboardInterrupt"
		pool.close()
		pool.terminate()
		pool.join()
		exit()

	lis_vcf_files = []
	for f in lis_bam_files:
		tmp_f = f.split("/")[-1]
		tmp_f = tmp_f[:-3]+"gvcf"
		lis_vcf_files.append(vcf_file_path+tmp_f)

	return lis_vcf_files

def get_gt(f1):
	fr = open(f1.file_path,"r")
	for line in fr:
		if line.startswith("#"):
			continue
		dp_idx = -1
		gt_idx = -1
		lis = line.strip().split("\t")
		lis_fm = lis[8].strip().split(":")
		for i in range(0,len(lis_fm)):
			if lis_fm[i] == "DP":
				dp_idx = i
			elif lis_fm[i] == "GT":
				gt_idx = i
		lis_info = lis[9].strip().split(":")
		if int(lis_info[dp_idx]) < 5:
			continue
		lis_gt = lis_info[gt_idx].split('/')
		try:
			gt = (int(lis_gt[0])+int(lis_gt[1]))*0.5
		except:
			lis_gt = lis_info[gt_idx].split('|')
			gt = (int(lis_gt[0])+int(lis_gt[1]))*0.5
		f1.dic_gt[lis[0]][lis[1]] = gt
	fr.close()

def cal_cor_each(f1,f2,cor_arr,x):
	count_all=0
	count_match=0
	count_diff_gt = 0
	lis_diff_snp = []
	for ch in f1.dic_gt.keys():
		for pos in f1.dic_gt[ch].keys():
			gt1 = f1.dic_gt[ch][pos]
			try:
				gt2 = f2.dic_gt[ch][pos]
			except:
				continue
			count_all += 1
			if gt1 == gt2:
				count_match += 1
			else:
				count_diff_gt += 1
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

def get_sw_pairs(lis_files,smp_pairs):
	separators = r'_|-|\.'
	dic_split_samples = {}
	dic_sample_scores = {}
	dic_topscores = {}
	lis_factorsize = []
	for v in lis_files:
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
	for f1 in lis_files:
		try:
			if smp_pairs[f1] == []:
				dic_un_p[f1] = get_max(dic_sample_scores[f1])
			else:
				dic_sw[f1] = []
				lis_m_f = get_max(dic_sample_scores[f1])
				for f2 in smp_pairs[f1]:
					if f2 not in lis_m_f:
						dic_sw[f1].append(f2)
				for f in lis_m_f:
					if f not in smp_pairs[f1]:
						dic_sw[f1].append(f)
		except:
			lis_m_f = get_max(dic_sample_scores[f1])
			for f in smp_pairs.keys():
				if f1 in smp_pairs[f]:
					if f not in lis_m_f:
						dic_sw[f1] = []
						for f_m in lis_m_f:
							dic_sw[f1].append(f_m)
	lis_sw_keys = dic_sw.keys()
	lis_sw_keys.sort()
	for f1 in lis_sw_keys:
		dic_sw[f1].sort()
		for f2 in dic_sw[f1]:
			try:
				if (f1 != f2) & (f1 in dic_sw[f2]):
					dic_sw[f2].remove(f1)
			except:
				pass
	return dic_sw, dic_un_p

def get_sw_pairs_ans(lis_files, smp_pairs, lis_ans):
	lis_f1 = smp_pairs.keys()
	lis_f1.sort()
	dic_sw = {}
	dic_un_p = {}
	for f1 in lis_f1:
		dic_sw[f1] = []
		dic_un_p[f1] = []
		for g in lis_ans:
			if f1 in g:
				g_f = copy.deepcopy(g)
				g_f.remove(f1)
				smp_pairs[f1].sort()
				if smp_pairs[f1] == []:
					dic_un_p[f1] = g_f
				if g_f != smp_pairs[f1]:
					if smp_pairs[f1] == []:
						dic_un_p[f1] = g_f
					else:
						if g_f[0] not in smp_pairs[f1]:
							dic_sw[f1].append(g_f[0])
						for f2 in smp_pairs[f1]:
							if f2 != g_f[0]:
								dic_sw[f1].append(f2)
								for g2 in lis_ans:
									if f2 in g2:
										g2_f = copy.deepcopy(g2)
										g2_f.remove(f2)
										try:
											dic_sw[f2].append(g2_f[0])
										except:
											dic_sw[f2]=[g2_f[0]]
	lis_sw_keys = dic_sw.keys()
	lis_sw_keys.sort()
	for f1 in lis_sw_keys:
		dic_sw[f1].sort()
		for f2 in dic_sw[f1]:
			if (f1 != f2) & (f1 in dic_sw[f2]):
				dic_sw[f2].remove(f1)
	return dic_sw, dic_un_p

def make_result_file_no_file_name_info(cor_matrix,smp_pairs,lis_files,OutputDIR,lis_paired_files):
	print "           Skip making 'Mismatched_samples.txt' file"
	fw_m_m = open(OutputDIR+"Matched_samples.txt","w")
	fw_m_m.write("# Matched pair by genotype\n")
	for i in range(0,len(lis_paired_files)):
		f1 = lis_paired_files[i]
		if smp_pairs[f1]  == []:
			continue
		for f2 in smp_pairs[f1]:
			flag_less_informative = False
			fw_m_m.write(f1+"\t"+f2+"\t")
			score = cor_matrix[lis_files.index(f1)][lis_files.index(f2)]
			if score < 0.8:
				flag_less_informative = True
			fw_m_m.write(str(score)+"\tMathced\n")
			if flag_less_informative:
				fw_m_m.write("-> *This pair scores under 0.8 which is less informative. The 'less informative score' dosen't mean that the pair is not matched but may have some problem of purity or copy number variation etc. in the sample.\n")
	fw_m_m.close()

def make_result_file(cor_matrix,smp_pairs,lis_files,OutputDIR,lis_ans):
	return_v = 1
	fw_a_m = open(OutputDIR+"Total_result.txt","w")
	lis_paired_files = smp_pairs.keys()
	lis_paired_files.sort()
	len_v = len(lis_files)
	lis_sw = []
	lis_up = []
	if ( lis_ans == [] ) & (len(lis_files) < 6) :
		print "## WARNING : The number of files is not enough to pair by file names."
		make_result_file_no_file_name_info(cor_matrix,smp_pairs,lis_files,OutputDIR,lis_paired_files)
		for i in range(0,len_v-1):
			for j in range(i+1, len_v):
				m_um="Unmatched"
				if cor_matrix[i][j] > 0.7:
					m_um = "Matched"
				fw_a_m.write(lis_files[i]+"\t"+lis_files[j]+"\t"+str(cor_matrix[i][j])+"\t"+m_um+"\n")
		fw_a_m.close()
		mk_html_no_mismatched(OutputDIR)
		return_v = 0
	else:
		count_m = 0
		count_s = 0
		count_u = 0
		dic_sw = {}
		dic_un_p = {}
		if lis_ans == []:
			dic_sw, dic_un_p = get_sw_pairs(lis_files,smp_pairs)
			if ( dic_sw == None ) & ( dic_un_p == None):
				print "## WARNING : The file names don't have detectable common regulation."
				make_result_file_no_file_name_info(cor_matrix,smp_pairs,lis_files,OutputDIR,lis_ans,lis_paired_files)
				mk_html_no_mismatched(OutputDIR)
				return_v = 0
				for i in range(0,len_v-1):
					for j in range(i+1, len_v):
			 			m_um="Unmatched"
						if cor_matrix[i][j] > 0.7:
							m_um = "Matched"
						fw_a_m.write(lis_files[i]+"\t"+lis_files[j]+"\t"+str(cor_matrix[i][j])+"\t"+m_um+"\n")
				fw_a_m.close()
				return return_v
			else:
				print "Detected pairs by file names."
		else:
			dic_sw, dic_un_p = get_sw_pairs_ans(lis_files, smp_pairs, lis_ans)
		fw_m_m = open(OutputDIR+"Matched_samples.txt","w")
		fw_s_m = open(OutputDIR+"Mismatched_samples.txt","w")
		perfect_m = True
		for f1 in lis_paired_files:
			if len(smp_pairs[f1]) > 0:
				for sw in dic_sw[f1]:
					perfect_m = False
					break
			else:
				for un in dic_un_p[f1]:
					perfect_m = False
					break
		if perfect_m:
			fw_s_m.write("No mismatched samples.")
		else:
			return_v = 2
		fw_m_m.write("# Matched pair by both genotype and name\n")
		lis_m = []
		for i in range(0,len(lis_paired_files)):
			f1 = lis_paired_files[i]
			if smp_pairs[f1]  == []:
				continue
			if len(smp_pairs[f1]) > 1:
				pass
			for f2 in smp_pairs[f1]:
				if f2 not in dic_sw[f1]:
					flag_less_informative = False
					fw_m_m.write(f1+"\t"+f2+"\t")
					score = cor_matrix[lis_files.index(f1)][lis_files.index(f2)]
					if score < 0.8:
						flag_less_informative = True
					fw_m_m.write(str(score)+"\tMatched\n")
					count_m += 1
					if flag_less_informative:
						fw_m_m.write("-> *This pair scores under 0.8 which is less informative. The 'less informative score' dosen't mean that the pair is not matched but may have some problem of purity or copy number variation etc. in the sample.\n")
					lis_m.append([f1,f2,round(score,2),"Matched"])
		fw_m_m.close()
		if not perfect_m:
			fw_s_m.write("# Matched samples only by genotype or file name but not by both\n")
			lis_sw_keys = dic_sw.keys()
			lis_sw_keys.sort()
			for f1 in lis_sw_keys:
				for f2 in dic_sw[f1]:
					flag_less_informative = False
					count_s +=  1
					fw_s_m.write(f1+"\t"+f2+"\t")
					score = cor_matrix[lis_files.index(f1)][lis_files.index(f2)]
				 	m_um="Unmatched"
					if score > 0.7 :
						if score < 0.8:
							flag_less_informative = True
						fw_s_m.write(str(score)+"\tMatched\n")
						m_um="Matched"
					else:
						fw_s_m.write(str(score)+"\tUnmathced\n")
					if flag_less_informative:
						fw_s_m.write("-> *This pair scores under 0.8 which is less informative. The 'less informative score' dosen't mean that the pair is not matched but may have some problem of purity or copy number variation etc. in the sample.\n")
					lis_sw.append([f1,f2,round(score,2),m_um])
			fw_s_m.write("\n# List of samples are matched with nothing by genotype\n")
			for i in range(0,len(lis_paired_files)):
				f1 = lis_paired_files[i]
				if smp_pairs[f1]  == []:
					count_u += 1
					fw_s_m.write(f1+"\n")
					for f2 in dic_un_p[f1]:
						score = cor_matrix[lis_files.index(f1)][lis_files.index(f2)]
						fw_s_m.write("-> pair by file name with "+ f2 +" ( score : "+str(score) +" )\n")
						lis_up.append([f1,f2,round(score,2),"Unmatched"])
		fw_s_m.close()

		count_line = 0
		for i in range(0,len_v-1):
			for j in range(i+1, len_v):
				count_line += 1
				try:
					if lis_files[j] in dic_sw[lis_files[i]]:
						for k in range(0,len(lis_sw)):
							sw = lis_sw[k]
							if sw[0] == lis_files[i]:
								if sw[1] == lis_files[j]:
									lis_sw[k].append(str(count_line+1))
				except:
					pass
				try:
					if lis_files[j] in dic_un_p[lis_files[i]]:
						for k in range(0,len(lis_up)):
							up = lis_up[k]
							if up[0] == lis_files[i]:
								if up[1] == lis_files[j]:
									lis_up[k].append(str(count_line+1))
							if up[0] == lis_files[j]:
								if up[1] == lis_files[i]:
									lis_up[k].append(str(count_line+1))
				except:
					pass
				m_um="Unmatched"
				if cor_matrix[i][j] > 0.7:
					m_um = "Matched"
				fw_a_m.write(lis_files[i]+"\t"+lis_files[j]+"\t"+str(cor_matrix[i][j])+"\t"+m_um+"\n")
	
		mk_html_dic(OutputDIR,lis_m,lis_sw,lis_up)
		fw_a_m.close()
	return return_v

def mk_html_dic(OutputDIR,lis_m,lis_sw,lis_up):
	fw_r = open(OutputDIR + "BAMixChecker_Report.Rmd","w")
	fw_r.write("# Sample Mix-up analysis result by BAMixChecker\n")
	fw_r.write("```{r, echo=FALSE , results='hide', message=FALSE, warning=FALSE}\n")
	fw_r.write("if(!(require(ztable))){install.packages('ztable')}\n")
	fw_r.write("library('ztable')\n")
	fw_r.write("dataDir='{0}'\n".format(OutputDIR))
	#matched
	lis_f1 = []
	lis_f2 = []
	lis_score = []
	lis_m_um = []
	for i in range(0,len(lis_m)):
		lis_f1.append("'"+lis_m[i][0]+"'")
		lis_f2.append("'"+lis_m[i][1]+"'")
		lis_score.append("'"+str(lis_m[i][2])+"'")
		lis_m_um.append("'"+lis_m[i][3]+"'")
	f1s = ','.join(lis_f1)
	f2s = ','.join(lis_f2)
	scores = ','.join(lis_score)
	m_ums = ','.join(lis_m_um)
	fw_r.write("df.m <-data.frame('Sample1'=c({0}), 'Sample2'=c({1}),'Concordance rate'=c({2}), 'Conclusion'=c({3}))\n".format(f1s,f2s,scores,m_ums))
	fw_r.write("colnames(df.m) <- c('Sample1', 'Sample2','Concordance rate', 'Conclusion')\n")
	#swapped
	if lis_sw != []:
		lis_f1 = []
		lis_f2 = []
		lis_score = []
		lis_m_um = []
		lis_sw_c = []
		for i in range(0,len(lis_sw)):
			lis_f1.append("'"+lis_sw[i][0]+"'")
			lis_f2.append("'"+lis_sw[i][1]+"'")
			lis_score.append("'"+str(lis_sw[i][2])+"'")
			lis_m_um.append("'"+lis_sw[i][3]+"'")
			lis_sw_c.append(lis_sw[i][4])
		f1s = ','.join(lis_f1)
		f2s = ','.join(lis_f2)
		scores = ','.join(lis_score)
		m_ums = ','.join(lis_m_um)
		fw_r.write("df.sw <-data.frame('Sample1'=c({0}), 'Sample2'=c({1}),'Concordance rate'=c({2}), 'Conclusion'=c({3}))\n".format(f1s,f2s,scores,m_ums))
		fw_r.write("colnames(df.sw) <- c('Sample1', 'Sample2','Concordance rate', 'Conclusion')\n")
	#unpaired
	if lis_up != []:
		lis_f1 = []
		lis_f2 = []
		lis_score = []
		lis_m_um = []
		lis_up_c = []
		for i in range(0,len(lis_up)):
			lis_f1.append("'"+lis_up[i][0]+"'")
			lis_f2.append("'"+lis_up[i][1]+"'")
			lis_score.append("'"+str(lis_up[i][2])+"'")
			lis_m_um.append("'"+lis_up[i][3]+"'")
			lis_up_c.append(lis_up[i][4])
		lis_up_c = list(set(lis_up_c))
		lis_up_c.sort()
		f1s = ','.join(lis_f1)
		f2s = ','.join(lis_f2)
		scores = ','.join(lis_score)
		m_ums = ','.join(lis_m_um)
		fw_r.write("df.up <-data.frame('Orphan sample'=c({0}), 'Best match by file name'=c({1}),'Concordance rate'=c({2}), 'Conclusion'=c({3}))\n".format(f1s,f2s,scores,m_ums))
		fw_r.write("colnames(df.up) <- c('Orphan sample', 'Best match by file name','Concordance rate', 'Conclusion')\n")
	fw_r.write("```\n## Mismatched samples\n")
	fw_r.write("#### 1. Swapped samples\n")
	fw_r.write("###### - Matched samples only by *genotype* or *file name* but not by both\n")
	if lis_sw != []:
		fw_r.write("```{r , results='asis', echo=FALSE}\n")
		fw_r.write("z = ztable(df.sw,align='cccc',include.rownames=FALSE)\n")
		fw_r.write("z <-addRowColor(z,c(1),'pink')\nprint (z, type = 'html')\n```\n")
	else:
		fw_r.write("##### *No swapped samples*\n")
	fw_r.write("#### 2. Orphan samples\n")
	fw_r.write("###### - Samples matched with nothing by genotype\n")
	if lis_up != []:
		fw_r.write("```{r , results='asis', echo=FALSE}\n")
		fw_r.write("z = ztable(df.up,align='cccc',include.rownames=FALSE)\n")
		fw_r.write("z <-addRowColor(z,c(1),'light green')\nprint (z, type = 'html')\n```\n")
	else:
		fw_r.write("##### *No unpaired samples*\n")
	fw_r.write("## Matched samples\n")
	fw_r.write("###### - The matched samples by *the file name* and *the genotype*\n")
	if lis_m != []:
		fw_r.write("```{r , results='asis', echo=FALSE}\n")
		fw_r.write("z = ztable(df.m,align='cccc',include.rownames=FALSE)\nprint (z, type = 'html')\n```\n")
	else:
		fw_r.write("##### *No matched samples*\n")
		
	fw_r.write("## Total result\n")
	fw_r.write("```{r , results='asis', echo=FALSE}\n")
	fw_r.write("df.total = read.delim(paste0(dataDir, 'Total_result.txt'), header=F)\n")
	fw_r.write("colnames(df.total) <- c('Sample1', 'Sample2','Concordance rate', 'Conclusion')\n")
	fw_r.write("z = ztable(df.total,align='cccc',include.rownames=FALSE)\n")
	if lis_sw != []:
		if lis_up != []:
			fw_r.write("z <-addRowColor(z,c({0}),'pink')\nz <-addRowColor(z,c({1}),'light green')\n".format(','.join(lis_sw_c),','.join(lis_up_c)))
		else:
			fw_r.write("z <-addRowColor(z,c({0}),'pink')\n".format(','.join(lis_sw_c)))
	else:
		if lis_up != []:
			fw_r.write("z <-addRowColor(z,c({0}),'light green')\n".format(','.join(lis_up_c)))
	fw_r.write("print (z, type = 'html')\n```\n")
	fw_r.close()
	cmm = "Rscript {0} {1} {2}".format(os.path.dirname(os.path.realpath(__file__))+"/r_script.r",OutputDIR + "BAMixChecker_Report.Rmd", OutputDIR)
	#os.system(cmm)
	prc= Popen(cmm, stdout=PIPE, shell=True, stderr=PIPE)
	stdoutput, stderr = prc.communicate()
	if prc.returncode != 0:
		print stderr

def	mk_html_no_mismatched(OutputDIR):
	fw_r = open(OutputDIR + "BAMixChecker_Report.Rmd","w")
	fw_r.write("# Sample Mix-up analysis result by BAMixChecker\n")
	fw_r.write("```{r, echo=FALSE , results='hide', message=FALSE, warning=FALSE}\n")
	fw_r.write("if(!(require(ztable))){install.packages('ztable')}\n")
	fw_r.write("library('ztable')\n")
	fw_r.write("dataDir='{0}'\n".format(OutputDIR))
	fw_r.write("df.total = read.delim(paste0(dataDir, 'Total_result.txt'), header=F)\n")
	fw_r.write("colnames(df.total) <- c('Sample1', 'Sample2','Concordance rate', 'Conclusion')\n")
	fw_r.write("```\n## Total result\n")
	fw_r.write("```{r , results='asis', echo=FALSE}\n")
	fw_r.write("z = ztable(df.total,align='llcl',include.rownames=FALSE)\nprint (z, type = 'html')\n```\n")
	fw_r.close()
	cmm = "Rscript {0} {1} {2}".format(os.path.dirname(os.path.realpath(__file__))+"/r_script.r",OutputDIR + "BAMixChecker_Report.Rmd", OutputDIR)
	prc= Popen(cmm, stdout=PIPE, shell=True, stderr=PIPE)
	stdoutput, stderr = prc.communicate()
	if prc.returncode != 0:
		print stderr


if __name__ == "__main__":
	start_t = time.time()
	parser = argparse.ArgumentParser(prog="BAMixChecker", description="Sample mix-up checker to detect sample mismatch with pairs of BAM file in a cohort for WGS/WES/RNA-seq and targeted sequencing.")
	parser.add_argument('-d','--DIR', default="", help="Directory path of the .BAM files")
	parser.add_argument('-l', '--List', default="", help="A file with the list of .BAM files")
	parser.add_argument('-r', '--Ref', default="",required=True, help="Reference file")
	parser.add_argument('-o','--OutputDIR', default="", help="Output directory path")
	parser.add_argument('-b','--BEDfile', default="", help="BED file for Targeted sequencing data")
	parser.add_argument('-v', '--RefVer', default="hg38", choices=['hg38','hg19'], help="Version of reference : 'hg19' or 'hg38'. Default = 'hg38'")
	parser.add_argument('-p', '--MaxProcess', default="4", help="The number of max process. Default = 4")
	parser.add_argument('--FullPATH', action='store_true',help="Use to report with the full path of file. BAMixChecker resports with the only file name as a default.")
	parser.add_argument('--RemoveVCF',action='store_true', help="Use to remove called germline VCF after running.")

	# get the tool path
	bedtools_path = ''
	HC_path = ""
	fr_config = open(os.path.dirname(os.path.realpath(__file__))+'/BAMixChecker.config','r')
	for line in fr_config:
		if line.startswith('BEDTOOLS'):
			bedtools_path = line.split('=')[1].strip()
		if line.startswith('GATK'):
			HC_path = line.split('=')[1].strip()
	fr_config.close()
	if bedtools_path == '':
		print "## ERROR: bedtools path is not set in 'BAMixChecker.config' file."
		exit()
	elif HC_path == "":
		print "## ERROR: GATK path is not set in 'BAMixChecker.config' file."
		exit()

	# get the arguments
	dir_path = ''
	out_path = ''
	args = parser.parse_args()
	if args.DIR == '':
		if args.List == '':
				print "## ERROR: There is no information about input files. Use -d or -l option for the input file information."
				exit()
	else:
		if args.List != '':
			print dir_path
			print "## ERROR: Option -d and -l are exclusive. Try with one of the options.\n##\t Check 'https://github.com/heinc1010/BAMixChecker' for more information about the options of BAMixChecker."
			exit()
		dir_path = abspath(args.DIR)+'/'
	if args.OutputDIR == '':
		out_path = abspath('.')+'/BAMixChecker/'
	else:
		out_path = abspath(args.OutputDIR)+'/BAMixChecker/'
	cmm = "mkdir -p {0}".format(out_path)
	os.system(cmm)

	if args.Ref == "":
		print "## ERROR: Reference file is necessary. Use -r option."
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
	if args.RefVer in [ "hg38","hg19" ]:
		if args.BEDfile != '':
			print "Run for targeted sequecing data"
			bed_file = make_bed("{0}gnomad_{1}_AF{2}_AF{3}_All.bed".format(bed_file_path,args.RefVer,45,35), args.BEDfile , out_path, bedtools_path)
			for AF in range(45,5,-5):
				for AF_all in range(AF,-1,-10):
					if AF_all >= AF:
						continue
					AF_all = int(AF_all/10)*10
					if AF_all != 0:
						bed_file = make_bed("{0}gnomad_{1}_AF{2}_AF{3}_All.bed".format(bed_file_path,args.RefVer,AF,AF_all), args.BEDfile , out_path, bedtools_path)
					else:
						bed_file = make_bed("{0}gnomad_{1}_AF{2}.bed".format(bed_file_path,args.RefVer,AF), args.BEDfile , out_path, bedtools_path)
					if bed_file != None:
						break
				if bed_file != None:
					break
		else:
			bed_file = "{0}gnomad_{1}_AF45_AF35_All.bed".format(bed_file_path,args.RefVer)
	else:
		print "## ERROR: Option -v should be 'hg19' or 'hg38'."
		exit()

	lis_bam_files,lis_ans = get_file_list(dir_path, args.List)

	if len(lis_bam_files) == 0:
		print "## ERROR: No .bam file in the list or directory."
		exit()

	# call the variants
	lis_vcf_files = run_HC(lis_bam_files,out_path,args.Ref,bed_file,args.MaxProcess,HC_path)

	# calculate the concordance
	cor_matrix = cal_cor(lis_vcf_files, args.MaxProcess)

#	# pair based on the genotype concordance 
	if args.FullPATH:
		# pair based on the genotype concordance 
		smp_pairs = pairing(cor_matrix,lis_bam_files)
		# determine the matched or mismatched pair based on file names as well as the genotype concordance
		result = make_result_file(cor_matrix,smp_pairs,lis_bam_files,out_path,lis_ans)
	else:
		lis_bam_files_sp = [ f.split("/")[-1] for f in lis_bam_files ]
		lis_ans_sp = []
		for ans in lis_ans:
			lis_ans_sp.append([f.split("/")[-1] for f in ans])
		# pair based on the genotype concordance 
		smp_pairs = pairing(cor_matrix,lis_bam_files_sp)
		# determine the matched or mismatched pair based on file names as well as the genotype concordance
		result = make_result_file(cor_matrix,smp_pairs,lis_bam_files_sp,out_path,lis_ans_sp)

	if result == 1:
		print "Perfect match."
	elif result == 2:
		print "Swapped file exist. Check 'BAMixChecker_Report.html' or 'Mismatched_samples.txt' file."

	if args.RemoveVCF:
		cmm = "rm -r {0}HaplotypeCaller/".format(out_path)
		os.system(cmm)

	end_t = time.time() -  start_t
	print "Running time: "+str(round(end_t/60,2))+" min" 
