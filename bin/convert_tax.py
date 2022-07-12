#!/usr/bin/env python

import re
import sys

tax_file = sys.argv[1]
out = sys.argv[2]

tax = open(tax_file,"r")
sure_tax_out = open(f"sure_{out}","w")
notSure_tax_out = open(f"notSure_{out}","w")

sure_tax_out.write("#OTU\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIE\n")
notSure_tax_out.write("#OTU\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIE\n")

for line in tax :
	otu=line.strip().split('\t')[0]
	try :
		sure_taxon_list=line.strip().split('\t')[3].split(',')
	except IndexError : #taxonomy not sure (under the threshold)
		sure_taxon_list=["d:NA","p:NA","c:NA","o:NA","f:NA","g:NA","s:NA"]
		print(line)
	
	try : 
		ns_taxon_list=line.strip().split('\t')[1].split(',')
	except IndexError :
		ns_taxon_list=["d:NA(1.00)","p:NA(1.00)","c:NA(1.00)","o:NA(1.00)","f:NA(1.00)","g:NA(1.00)","s:NA(1.00)"]

	sure_strip_taxon_l = []
	for s_t in sure_taxon_list :
		sure_taxon = re.search(':(.*)',s_t)
		if sure_taxon.group(1)[0] == '"' :
			sure_taxon = sure_taxon.group(1).strip("\"")
		else :
			sure_taxon = sure_taxon.group(1)
		sure_strip_taxon_l.append(sure_taxon)

	ns_strip_taxon_l = []
	for ns_t in ns_taxon_list :
		ns_taxon = re.search(':(.*)\(',ns_t)
		if ns_taxon.group(1)[0] == '"' :
			ns_taxon = ns_taxon.group(1).strip("\"")
		else :
			ns_taxon = ns_taxon.group(1)
		ns_strip_taxon_l.append(ns_taxon)


	sure_tax_out.write(f"{otu}\t")
	notSure_tax_out.write(f"{otu}\t")

	for sure_taxon in sure_strip_taxon_l :
		sure_tax_out.write(f"{sure_taxon}\t")
	sure_tax_out.write("\n")

	for ns_taxon in ns_strip_taxon_l :
		notSure_tax_out.write(f"{ns_taxon}\t")
	notSure_tax_out.write("\n")
