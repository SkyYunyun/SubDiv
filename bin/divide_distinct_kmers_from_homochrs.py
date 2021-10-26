#!/usr/bin/python3
import sys
trantab=str.maketrans('ABCDEFGHIGKLMNOPQRSTUVWXYZabcdefghigklmnopqrstuvwxyz', 'TBGDEFCHIGKLMNOPQRSAUVWXYZtbgdefghicklmnopqrsauvwxyz')
def DNA_complements(seq):
	string=seq.translate(trantab)
	return string[::-1]
def read_seq_as_dict(fastafile):
	seqdict={}
	with open(fastafile,'r') as fi:
		for line in fi:
			line=line.strip()
			if line.startswith('>'):
				ID=line.split()[0].replace('>','')
				seqdict[ID]=[]
			else:
				seqdict[ID].append(line)
		for ID in seqdict:
			seqdict[ID]=''.join(seqdict[ID])
			seqdict[ID]=seqdict[ID].upper()
	return seqdict
def homo_kmer_pair_align(homo_groups,sequence_dict,kmer_result,kmer_len=13,min_times=2):
	kmerDictList=[{} for i in range(len(homo_groups))]
	for ID_index in range(len(homo_groups)):
		sequence=sequence_dict[homo_groups[ID_index]]
		i=0;seq_len=len(sequence)
		while i<=seq_len-kmer_len:
			kmer=sequence[i:i+kmer_len]
			i+=1
			if 'N' in kmer:continue
			com_rev_kmer=DNA_complements(kmer)
			if com_rev_kmer > kmer: kmer=com_rev_kmer
			if kmer not in kmerDictList[ID_index]:
				kmerDictList[ID_index][kmer]=0
			kmerDictList[ID_index][kmer]+=1
	all_kmers_set=set()
	for kmerDict in kmerDictList:
		for kmer in kmerDict:
			all_kmers_set.add(kmer)
	header=['Kmer']+homo_groups
	headerline='\t'.join(header)
	fo=open(kmer_result,'w')
	print(headerline,file=fo)
	for kmer in all_kmers_set:
		homo_counts=[]
		for kmerDict in kmerDictList:
			if kmer in kmerDict:
				homo_counts.append(kmerDict[kmer])
			else:
				homo_counts.append(0)
		if max(homo_counts)-min(homo_counts)>=min_times*min(homo_counts) and max(homo_counts)>=min_times:
			linelist=[kmer]+[str(count) for count in homo_counts]
			line='\t'.join(linelist)
			print(line,file=fo)
def mutiple_cpu_align_and_combine(sequence_dict,homo_chr_file,cpu_num=10,kmer_len=13,min_times=2):
	import multiprocessing
	from multiprocessing import Pool
	process_pool=Pool(cpu_num);outfiles=[];
	with open(homo_chr_file,'r') as fi:
		n=0
		for line in fi:
			line=line.strip()
			homo_groups=line.split()
			outfilename='g'+str(n+1)+'_kmer'
			outfiles.append(outfilename)
			process_pool.apply_async(homo_kmer_pair_align,args=(homo_groups,sequence_dict,outfilename,kmer_len,min_times))
			sys.stdout.write("put task to process_pool!\n")
			sys.stdout.flush()
			n+=1
	process_pool.close()
	process_pool.join()
	sys.stdout.write("loading kmers complete!\n")
	sys.stdout.flush()
	kmerDictList=[{}for i in range(len(outfiles))]
	kmerset=set()
	for file_index in range(len(outfiles)):
		with open(outfiles[file_index],'r') as fi:
			for line in fi:
				line=line.strip()
				elements=line.split()
				kmerDictList[file_index][elements[0]]=elements[1:]
				kmerset.add(elements[0])
	signtrue='T'*len(outfiles)
	finalfile=open('distictive_kmer_and_counts','w')
	kmerlist=list(kmerset)
	kmerlist.remove('Kmer')
	kmerlist=['Kmer']+kmerlist
	for kmer in kmerlist:
		sign_list=[]
		lines=[]
		for kmerDict in kmerDictList:
			if kmer in kmerDict:
				sign_list.append('T')
				lines.extend(kmerDict[kmer])
			else:
				sign_list.append('F')
		lines=[kmer]+lines
		if ''.join(sign_list)==signtrue:
			print('\t'.join(lines),file=finalfile)
if __name__=='__main__':
	genome_fasta=sys.argv[1]
	homolog_groups=sys.argv[2]
	cpu_num=int(sys.argv[3])
	kmer_len=int(sys.argv[4])
	min_times=int(sys.argv[5])
	seqdict=read_seq_as_dict(genome_fasta)	
	mutiple_cpu_align_and_combine(seqdict,homolog_groups,cpu_num=cpu_num,kmer_len=kmer_len,min_times=min_times)
