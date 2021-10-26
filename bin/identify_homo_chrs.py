#!/usr/bin/env python3

total=list(range(32))
def average_divide_list(total_list,n_parts):
	n_part_list=[]
	total_len=len(total_list)
	step=int(total_len/n_parts)+1
	for i in range(0,total_len,step):
		n_part_list.extend([total_list[i:i+step]])
	return n_part_list
def read_fasta_as_dict(fastafile):
	fasta_dict={}
	with open(fastafile,'r') as fi:
		for line in fi:
			line=line.strip()
			if line.startswith('>'):
				ID=line.split()[0].replace('>','')
				fasta_dict[ID]=[]
			else:
				fasta_dict[ID].append(line.upper())
		for ID in fasta_dict:
			fasta_dict[ID]=''.join(fasta_dict[ID])
	return fasta_dict
def sliding_fasta(fasta_dict,step=100,window=100,outfile='split'):
	with open(outfile,'w') as fo:
		for ID in fasta_dict:
			sequence=fasta_dict[ID]
			seq_len=len(sequence)
			for i in range(0,seq_len,step):
				subseq=sequence[i:i+window]
				if len(subseq)==window:
					new_ID='>'+str(ID)+':'+str(i+1)+'-'+str(i+window)
					print(new_ID,file=fo)
					print(subseq,file=fo)
def write_fas(fasta_dict,target_list,filename):
	fo=open(filename,'w')
	for ID in target_list:
		fo.write('>%s\n%s\n' % (ID,fasta_dict[ID]))
def divide_fasta_nparts(fasta_dict,newfileprefix='split_',cpunumber=1):
	import multiprocessing;split_name_list=[]
	multiprocessing_pool=multiprocessing.Pool(cpunumber)
	name_list=list(fasta_dict.keys())
	#print(name_list)
	n_part_list=average_divide_list(name_list,cpunumber)
	#print(n_part_list[0:2])
	for i in range(len(n_part_list)):
		file_num=i+1
		file_name=newfileprefix+str(file_num)
		split_name_list.append(file_name)
		multiprocessing_pool.apply_async(write_fas,(fasta_dict,n_part_list[i],file_name))
	multiprocessing_pool.close()
	multiprocessing_pool.join()
	return split_name_list
def run_blastn(dbname,split_files,cpunumber=1,threading=6,outfmt=6,evalue=1e-5,max_target_seqs=6,perc_identity=80,qcov_hsp_perc=95):
	import os
	outfileslist=[]
	if not os.path.exists(dbname+'.nsq'):
		cmd='makeblastdb -in '+dbname+' -out '+dbname+' -dbtype nucl'
		os.system(cmd)
	import multiprocessing
	cpu_pool=multiprocessing.Pool(cpunumber)
	for f in split_files:
		outfileslist.append(f+'.blastn.fmt6')
		cmd='blastn -db '+dbname+' -query '+f+' -out '+f+'.blastn.fmt6 -evalue '+str(evalue)+' -num_threads '+str(threading)+' -outfmt '+str(outfmt)+' -max_target_seqs '+str(max_target_seqs)+' -perc_identity '+str(perc_identity)+' -qcov_hsp_perc '+str(qcov_hsp_perc)
		cpu_pool.apply_async(func=os.system,args=(cmd,))
	cpu_pool.close()
	cpu_pool.join()
	return outfileslist
def select_subbesthit(blastfile,outfile):
	with open(blastfile,'r') as fi:
		fo=open(outfile,'w')
		subbesthitlist=[]
		ID='';score=0
		for line in fi:
			line=line.strip()
			elements=line.split()
			if float(elements[2])==100:continue
			currentscore=float(elements[2])*(int(elements[7])-int(elements[6])+1)
			if elements[0]!=ID:
				ID=elements[0]
				subbesthitlist.append(line)
				score=currentscore
			else:
				if currentscore>score:
					subbesthitlist.pop()
					subbesthitlist.append(line)
					score=currentscore
		for e in subbesthitlist:
			fo.write('%s\n' % e)
		fo.close()
def mutipcocess_combine_blast_res(outfileslist,outputfile,cpunumber=1):
	subbesthitfilelist=[]
	import multiprocessing
	cpu_pool=multiprocessing.Pool(cpunumber)
	for f in outfileslist:
		fout=f+'.subbesthit'
		subbesthitfilelist.append(fout)
		cpu_pool.apply_async(func=select_subbesthit,args=(f,fout))
	cpu_pool.close()
	cpu_pool.join()
	fo=open(outputfile,'w')
	for f in subbesthitfilelist:
		for line in open(f,'r'):
			fo.write('%s' % line)
	fo.close()
def counts_subbesthit(subbestfile,countmatrixfile,pairsfile):
	matrixdict={}
	idset=set()
	with open(subbestfile,'r') as fi:
		for line in fi:
			line=line.strip()
			elements=line.split()
			ID=elements[0].split(':')[0]
			idset.add(elements[1])
			if elements[1] not in matrixdict:matrixdict[elements[1]]={}
			if ID not in matrixdict[elements[1]]:matrixdict[elements[1]][ID]=0
			matrixdict[elements[1]][ID]+=1
			if elements[1]==ID:matrixdict[elements[1]][ID]=0
	idlist=list(idset)
	with open(countmatrixfile,'w') as fo:
		foo=open(pairsfile,'w')
		header='\t'.join(['ID']+list(idset))
		fo.write('%s\n' % header)
		idset=set()
		for IDf in idlist:
			linelist=[]
			valuelist=[]
			for IDs in idlist:
				if IDs in matrixdict[IDf]:
					linelist.append(str(matrixdict[IDf][IDs]))
					valuelist.append(matrixdict[IDf][IDs])
				else:
					linelist.append('0')
					valuelist.append(0)
			maxpossible=max(valuelist)
			maxindex=idlist[valuelist.index(maxpossible)]
			tmplist=[]
			if IDf > maxindex:
				tmplist.extend([IDf,maxindex])
			else:
				tmplist.extend([maxindex,IDf])
			if ''.join(tmplist) not in idset:
				foo.write('%s\t%s\n' % (IDf,maxindex))
				idset.add(''.join(tmplist))
			line='\t'.join([IDf]+linelist)
			fo.write('%s\n' % line)
if __name__=='__main__':
	import sys
	fastafile=sys.argv[1]
	'''
	step1
	'''
	fasta_dict=read_fasta_as_dict(fastafile)
	step=int(sys.argv[2])
	window=int(sys.argv[3])
	cpunumber=int(sys.argv[4])
	outfile=fastafile.split('/')[-1]+'.split_step_'+str(step)+'_window_'+str(window)
	sliding_fasta(fasta_dict,step=step,window=window,outfile=outfile)
	split_fasta_dict=read_fasta_as_dict(outfile)
	newfileprefix=outfile+'.split_'
	split_name_list=divide_fasta_nparts(split_fasta_dict,newfileprefix=newfileprefix,cpunumber=cpunumber)
	blastnthreading=6
	blastnres_fileslist=run_blastn(fastafile,split_name_list,cpunumber=cpunumber,threading=blastnthreading,outfmt=6,evalue=1e-5,max_target_seqs=6,perc_identity=80,qcov_hsp_perc=90)
	mutipcocess_combine_blast_res(blastnres_fileslist,'subblasthit.res',cpunumber=cpunumber)
	counts_subbesthit('subblasthit.res','subblasthitcountmatrix.res','pairs_res_file')
