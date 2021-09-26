//
// Alignment with STAR
//

params.align_options          = [:]
params.samtools_sort_options  = [:]
params.samtools_index_options = [:]
params.samtools_stats_options = [:]

include { STAR_ALIGN        } from '../../modules/nf-core/modules/star/align/main' addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS } from '../nf-core/bam_sort_samtools'                           addParams( sort_options: params.samtools_sort_options, index_options: params.samtools_index_options, stats_options: params.samtools_stats_options )
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/modules/samtools/index/main' addParams( options: params.samtools_index_options )
include { BAM_STATS_SAMTOOLS } from '../nf-core/bam_stats_samtools'                              addParams( options: params.samtools_stats_options )
workflow ALIGN_STAR {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/star/index/
    gtf   // channel: /path/to/genome.gtf

    main:

    //
    // Map reads with STAR
    //
    STAR_ALIGN ( reads, index, gtf )

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    // Default output of STAR was set as sortedByCoord.out.bam, then no need to sort

    SAMTOOLS_INDEX     ( STAR_ALIGN.out.bam )
    STAR_ALIGN.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS ( ch_bam_bai )       
    
    emit:
    orig_bam         = STAR_ALIGN.out.bam             // channel: [ val(meta), bam            ]
    log_final        = STAR_ALIGN.out.log_final       // channel: [ val(meta), log_final      ]
    log_out          = STAR_ALIGN.out.log_out         // channel: [ val(meta), log_out        ]
    log_progress     = STAR_ALIGN.out.log_progress    // channel: [ val(meta), log_progress   ]
    bam_sorted       = STAR_ALIGN.out.bam_sorted      // channel: [ val(meta), bam_sorted     ]
    bam_transcript   = STAR_ALIGN.out.bam_transcript  // channel: [ val(meta), bam_transcript ]
    fastq            = STAR_ALIGN.out.fastq           // channel: [ val(meta), fastq          ]
    tab              = STAR_ALIGN.out.tab             // channel: [ val(meta), tab            ]
    star_version     = STAR_ALIGN.out.version         // path: *.version.txt

    bam              = STAR_ALIGN.out.bam             // channel: [ val(meta), [ bam ] ]
    bai              = SAMTOOLS_INDEX.out.bai         // channel: [ val(meta), [ bai ] ]
    csi              = SAMTOOLS_INDEX.out.csi         // channel: [ val(meta), [ csi ] ]
    stats            = BAM_STATS_SAMTOOLS.out.stats   // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]
    samtools_version = SAMTOOLS_INDEX.out.version        //    path: *.version.txt
}
