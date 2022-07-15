#!/usr/bin/env python

import re
import sys

# we input the .sintax file and name the output file
tax_file = sys.argv[1]
out = sys.argv[2]

tax = open(tax_file,"r") #open in read mode the tax file 
sure_tax_out = open(f"sure_{out}","w") #create an output file based on the name pass in parameters : it will contains taxonomy upper than the threshold put in parameters during the taxonomy
notSure_tax_out = open(f"notSure_{out}","w") # "" : it will contains all taxonomies, even if the score is under the threshold

sure_tax_out.write("#OTU\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIE\n") #write the headers
notSure_tax_out.write("#OTU\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIE\n")

for line in tax : # for each OTU/ASV
	otu=line.strip().split('\t')[0] # recover the name of the OTU/ASV
	#print(otu)
	try :
		sure_taxon_list=line.strip().split('\t')[3].split(',') #recover the list of taxons which have score upper the threshold
	except : # if existing no taxonomies upper than the threshold
		sure_taxon_list=["d:NA","p:NA","c:NA","o:NA","f:NA","g:NA","s:NA"] #put 'NA' value for each taxonomic ranks
	#print(sure_taxon_list)
	
	try : 
		ns_taxon_list=line.strip().split('\t')[1].split(',') # recover the list if taxon, even those under the threshold
		if ns_taxon_list==[''] :
			raise ValueError('No taxonomies found')
	except : # if not
		ns_taxon_list=["d:NA(1.00)","p:NA(1.00)","c:NA(1.00)","o:NA(1.00)","f:NA(1.00)","g:NA(1.00)","s:NA(1.00)"] #put 'NA' values for each ranks wih the maximum score
	#print(ns_taxon_list)		

	sure_strip_taxon_l = []
	for s_t in sure_taxon_list : # for each rank in all taxonomies of this OTU/ASV
		sure_taxon = re.search(':(.*)',s_t) # we recover only the exact taxon name (string after the :)
		if sure_taxon.group(1)[0] == '"' : # if the taxon start with a " 
			sure_taxon = sure_taxon.group(1).strip("\"") #remove it
		else :
			sure_taxon = sure_taxon.group(1)
		sure_strip_taxon_l.append(sure_taxon) #append this exact taxon name to the final list of sure taxons, for this OTU/ASV

	ns_strip_taxon_l = []
	for ns_t in ns_taxon_list : # for each taxo rank of this OTU/ASV
		ns_taxon = re.search(':(.*)\(',ns_t) # recover only the exact name of the taxon (string between : and ( )
		if ns_taxon.group(1)[0] == '"' : # if the taxon start with a "
			ns_taxon = ns_taxon.group(1).strip("\"") #remove it
		else :
			ns_taxon = ns_taxon.group(1)
		ns_strip_taxon_l.append(ns_taxon) #append this exact taxon name to the final list of taxons, for this OTU/ASV


	#write the name of the taxon in both output file
	sure_tax_out.write(f"{otu}\t") 
	notSure_tax_out.write(f"{otu}\t")

	# write each taxonomic rank in files 
	for sure_taxon in sure_strip_taxon_l :
		sure_tax_out.write(f"{sure_taxon}\t")
	sure_tax_out.write("\n")

	for ns_taxon in ns_strip_taxon_l :
		notSure_tax_out.write(f"{ns_taxon}\t")
	notSure_tax_out.write("\n")
