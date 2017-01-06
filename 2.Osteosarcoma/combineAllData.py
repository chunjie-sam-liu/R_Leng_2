#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name: combineAllData.py
#Author: C.J. Liu
#Mail: samliu@hust.edu.cn
#Created Time: Mon 07 Nov 2016 03:23:39 PM CST
################################################



import os,sys
import re
from glob import glob

filePath = '/extraspace/TARGET/OS/mRNA-seq/L3/expression/NCI-Meltzer'

def run():
	geneExpTxts = glob(filePath+"/*-01A-01R.gene.quantification.txt")
	# print(geneExpTxts)
	for txt in geneExpTxts:
		fileName = os.path.basename(txt)
		# print(fileName)
		name = re.sub(r'-01A-01R\.gene\.quantification\.txt','',fileName)
		# print(name)
		with open(txt,'r') as foo:
			pass

def main():
	run()

if __name__ == '__main__':
	main()

