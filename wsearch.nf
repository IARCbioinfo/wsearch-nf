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
params.output_folder	 = "Output_folder"
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
log.info "  Microbiome.nf 1.0 : Microbiome analysis "
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
    log.info "nextflow run iarcbioinfo/microbiome.nf --input_folder fastq/"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--input_folder        FOLDER               input folder with fastq files"
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

fastq_all_pairs = Channel.fromFilePairs(params.input_folder+'/*{L001_R1,L001_R2}_001.fastq')
fastq_all = Channel.fromPath(params.input_folder+'/*.fastq')

process cacheSilva {
  storeDir 'db/silva'

  output:
  file '*.fasta.gz' into silva
  file 'rdp_gold.fa' into chimera

  script:
  """
  makeSilvaRef.sh
  """
}

			
process run_usearch_summary {

     publishDir params.output_folder+"/1.reads_summary", mode: 'copy'

     input:
     file(fq) from fastq_all

     output:
     file '*info.txt' into summary
  
     shell:
     '''
     usearch -fastx_info !{fq} -output !{fq}.1a_fwd_fastq_info.txt > !{fq}.log
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
     usearch -fastq_mergepairs !{fq[0]} -reverse !{fq[1]} -fastqout !{sampleID}"_contig_!{params.experiment}.fastq" -fastq_maxdiffs !{params.maxdiff} -fastq_minovlen !{params.overlap} -report 2a_merging_seqs_report_!{params.experiment}.txt -tabbedout 2b_tabbedout_!{params.experiment}.txt -fastq_pctid !{params.identity} -fastq_minmergelen !{params.min_contig_length} -fastq_maxmergelen !{params.max_contig_length} > !{sampleID}.log
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
	usearch -fastq_filter !{fq} -fastaout !{ID}"_filtered_!{params.experiment}.fasta" -fastq_maxee !{params.max_ee} > !{ID}.log
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

        shell:
        '''
	usearch -fastx_uniques !{fa} -fastaout !{ID}"_uniq_!{params.experiment}.fasta" -sizeout
	'''
} 

process run_usearch_lowAbund {

        publishDir params.output_folder+"/6.low-abund_data", mode: 'copy'

        input:
        set val(ID), file(fa) from derep

        output:
        set val(ID), file('*low-abund*.fasta') into lowAbund

        shell:
        '''
	usearch -sortbysize !{fa} -fastaout !{ID}"_low-abund_!{params.experiment}.fasta" -maxsize 1
	'''
}

process run_usearch_exact {

        publishDir params.output_folder+"/7.Singleton-free", mode: 'copy'

        input:
        set val(ID), file(fa_la), file(fa_l) from lowAbund.join(labelled2)

        output:
        file('*SF*.fasta') into SF

        shell:
        '''
	usearch -search_exact !{fa_l} -db !{fa_la} -strand plus -notmatched !{ID}"_SF_!{params.experiment}.fasta"
	'''
}


process run_usearch_denoiseOtuAsv_taxonomy {

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

	shell:
        '''
        denoiseOtuAsv_taxonomy.sh !{fastaFiles} !{labelledFiles} !{chimera} !{silva} !{params.threads}
	'''
}

process run_usearch_Ranalysis {

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
