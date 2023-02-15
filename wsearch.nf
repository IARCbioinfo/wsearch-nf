#! /usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help          	 = null
params.input_folder  	 = "microbiome_analysis"
params.output_folder	 = "./"
params.maxdiff		 = "40"
params.experiment	 = "Experiment"
params.overlap		 = "50"
params.identity		 = "70"
params.min_contig_length = "50"
params.max_contig_length = "550"
params.max_ee		 = "1"
params.threads	         = "20"

log.info ""
log.info "----------------------------------------------------------------"
log.info "  Wsearch.nf 1.0 : Microbiome analysis with Usearch, Vsearch and Phyloseq"
log.info "  by Pierre BERTRAND : pierre.bertrand148@gmail.com"
log.info "----------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/wsearch.nf --input_folder fastq/ --output_folder output/ --experiment expirement_name"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--input_folder        FOLDER               Path to the folder containing .fastq files (forward and reverse) of each samples"
    log.info "--output_folder       FOLDER               Folder where the pipeline will perform and output"
    log.info "--experiment          STRING               Name of the experiment"
    log.info "" 
    log.info "Optionals arguments:"   
    log.info "--maxdiff             INT                  Maximum number of missmatch between the forward and the reverse read, in the overlap region. Can be decreased if reads have great quality [default=40]"
    log.info "--overlap             INT                  Minimum length (base pairs) of the overlap between the forward and the reverse reads. Couple of reads with a smaller overlap are removed [default=50]"
    log.info "--identity            INT                  Percentage of identity minimum between the overlap region of the forward and reverse reads. Can be increased with reads have good quality. [default=70]"
    log.info "--min_contig_length   INT                  Minimum length of the contig, after merge of the forward and reverse read [default=50]"
    log.info "--max_contig_length   INT                  Maximum length of the contig, after merge of the forward and reverse read [default=550]"
    log.info "--max_ee              NUM                  Maximum of expected errors allowed in each contigs. Expected errors are linked to the quality score from the fastq, and permit to estimate the number of errors in the contig. Can be increased if reads have a bad quality. https://drive5.com/usearch/manual8.0/expected_errors.html"
    log.info "--threads             INT                  Number of threads pipeline can use. Used only for steps need lot of computing power (Clustering, Denoising, Taxonomy, Analysis). [default=20]"
    log.info ""
    exit 0
}

/* Software information */
log.info ""
log.info "input_folder           	= ${params.input_folder}"
log.info "output_folder			= ${params.output_folder}"
log.info "maxdiff			= ${params.maxdiff}"
log.info "overlap                       = ${params.overlap}"
log.info "experiment	                = ${params.experiment}"
log.info "identity			= ${params.identity}"
log.info "min_contig_length		= ${params.min_contig_length}"
log.info "max_contig_length             = ${params.max_contig_length}"
log.info "max_ee			= ${params.max_ee}"
log.info "threads			= ${params.threads}"
log.info ""


fastq_allgz_pairs = Channel.fromFilePairs(params.input_folder+'/*{L001_R1,L001_R2}_001.fastq.gz')
//fastq_all = Channel.fromPath(params.input_folder+'/*.fastq.gz').view()


process cacheSilva {
  storeDir 'db/silva'

  output:
  file 'SILVA138_RESCRIPt.fasta' into silva
  file 'rdp_gold.fa' into chimera

  script:
  """
  makeSilvaRef.sh
  """
}

process decompress { 



     input:
     tuple val(sampleID), path(fq) from fastq_allgz_pairs

     output:
     file '*fastq' into fastq_all
     tuple val(sampleID), file('*_R{1,2}.fastq') into fastq_all_pairs

     shell:
     '''
     gzip -dc !{fq[0]} > !{sampleID}_R1.fastq
     gzip -dc !{fq[1]} > !{sampleID}_R2.fastq
     '''
}

			
process run_usearch_summary {


     publishDir params.output_folder+"/1.reads_summary", mode: 'copy'

     input:
     file(fq) from fastq_all.flatten()

     output:
     file 'SummaryReport*.log' into summary
  
     shell:
     '''
     usearch -fastx_info !{fq} -output !{fq}.1a_fwd_fastq_info.txt > SummaryReport_!{fq}.log
     '''
  }

process run_usearch_merge {


	publishDir params.output_folder+"/2.merge_reads", mode: 'copy'

     input:
     tuple val(sampleID), file(fq) from fastq_all_pairs

     output:
     set val(sampleID), file('*contig*.fastq') into contigs
     file('*.txt') into report_merge
     file('*.log') into log1

     shell:
     '''
     usearch -fastq_mergepairs !{fq[0]} -reverse !{fq[1]} -fastqout !{sampleID}_contig_!{params.experiment}.fastq -fastq_maxdiffs !{params.maxdiff} -fastq_minovlen !{params.overlap} -report 2a_merging_seqs_report_!{sampleID}_!{params.experiment}.txt -tabbedout 2b_tabbedout_!{sampleID}_!{params.experiment}.txt -fastq_pctid !{params.identity} -fastq_minmergelen !{params.min_contig_length} -fastq_maxmergelen !{params.max_contig_length} 2> MergeReport_!{sampleID}_!{params.experiment}.log
     '''
  }

process run_usearch_filter {

	
	publishDir params.output_folder+"/3.filtered_contigs", mode: 'copy'

	input: 
	set val(ID), file(fq) from contigs

	output:
	set val(ID), file('*filtered*.fasta') into filtered	
	file('*.log') into log2	

	shell:
	'''
	usearch -fastq_filter !{fq} -fastaout !{ID}"_filtered_!{params.experiment}.fasta" -fastq_maxee !{params.max_ee} 2> FilterReport_!{ID}_!{params.experiment}.log
	'''
}

process run_usearch_rename {


        publishDir params.output_folder+"/4.labeled_data", mode: 'copy'

        input:
	set val(ID), file(fa) from filtered

        output:
        set val(ID), file('*labelled*.fasta') into labelled
	file('*labelled*.fasta') into labelled3

        shell:
        '''
	sed -e "s/>/>barcodelabel=!{ID};/g" !{fa} > !{ID}"_labelled_!{params.experiment}.fasta"
	'''
}

process run_usearch_derep {


        publishDir params.output_folder+"/5.dereplicated_data", mode: 'copy'

        input:
        set val(ID), file(fa) from labelled

        output:
        set val(ID), file('*uniq*.fasta') into derep
	set val(ID), file(fa) into labelled2	
	file('*.log') into log5

        shell:
        '''
	usearch -fastx_uniques !{fa} -fastaout !{ID}"_uniq_!{params.experiment}.fasta" -sizeout > >(tee derep_!{ID}.log) 2> >(tee -a derep_!{ID}.log >&2)
	'''
} 

process run_usearch_lowAbund {


        publishDir params.output_folder+"/6.low-abund_data", mode: 'copy'

        input:
        set val(ID), file(fa) from derep

        output:
        set val(ID), file('*low-abund*.fasta') into lowAbund
	file('*.log') into log6

        shell:
        '''
	usearch -sortbysize !{fa} -fastaout !{ID}"_low-abund_!{params.experiment}.fasta" -maxsize 1  > >(tee lowAbund_!{ID}.log) 2> >(tee -a lowAbund_!{ID}.log >&2)
	'''
}

process run_usearch_exact {

        publishDir params.output_folder+"/7.Singleton-free", mode: 'copy'

        input:
        set val(ID), file(fa_la), file(fa_l) from lowAbund.join(labelled2)

        output:
        file('*SF*.fasta') into SF
	file('*.log') into log7

        shell:
        '''
	usearch -search_exact !{fa_l} -db !{fa_la} -strand plus -notmatched !{ID}"_SF_!{params.experiment}.fasta" > >(tee exact_!{ID}.log) 2> >(tee -a exact_!{ID}.log >&2)
	'''
}


process run_usearch_denoiseOtuAsv_taxonomy {

	label 'big_mem'

        publishDir params.output_folder+"/8.Denoising-Taxonomy/", mode: 'copy'
	cpus params.threads
		
        input:
        file(fastaFiles) from SF.collect()
        file(labelledFiles) from labelled3.collect()
        file(silva) from silva
        file(chimera) from chimera

        output:
        file('*.fasta') into fa
        file('8*') into denoising
	file('Taxonomy') into Taxonomy	
        file('8A.otu/uparse_otu_tab.txt') into otuTable
	file('Taxonomy/sure_sintax_newTax_otu.txt') into otuTax
	file('8B.asv/unoise_zotus_tab.txt') into asvTable
	file('Taxonomy/sure_sintax_newTax_asv.txt') into asvTax
	file('*.log') into log8

	shell:
        '''
	echo "\n ||| Pool files & Dereplication ||| \n"

	#pool all files in one
	cat !{fastaFiles} > all_SF.fasta
	cat !{labelledFiles} > all_labelled.fasta

        denoiseOtuAsv_taxonomy.sh all_SF.fasta all_labelled.fasta !{chimera} !{silva} !{params.threads}
	'''
}

process run_usearch_Ranalysis {

	label 'mid_mem'

        publishDir params.output_folder, mode: 'copy'
        cpus params.threads

        input:
	file(otuTable) from otuTable
	file(asvTable) from asvTable
	file(otuTax) from otuTax
	file(asvTax) from asvTax

        output:
        file('Analysis') into analyse

	shell:
	'''
	Usearch2phyloseq.R !{otuTable} !{otuTax} !{asvTable} !{asvTax}
	'''
}
