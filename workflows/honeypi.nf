include { DOWNLOAD_DB       } from '../modules/local/download_db/main'
include { TRIM_GALORE       } from '../modules/local/trim_galore/main'
include { DADA2             } from '../modules/local/dada2/main'
include { ITSX              } from '../modules/local/itsx/main'
include { CONSOLIDATE       } from '../modules/local/consolidate/main'
include { RDP_CLASSIFIER    } from '../modules/local/rdp_classifier/main'
include { FILTER_ASV_TABLE  } from '../modules/local/filter_asv_table/main'
include { MERGE_DUPLICATES  } from '../modules/local/merge_duplicates/main'

workflow HONEYPI {

    take:
    ch_reads   // [ meta, r1, r2 ]

    main:
    ch_versions = Channel.empty()

    // ── Database ─────────────────────────────────────────────────────────
    if (params.rdp_db_dir) {
        def db_url = (params.its_region == 'ITS1') ? params.rdp_db_its1_url : params.rdp_db_its2_url
        ch_db = Channel.value([file(params.rdp_db_dir), db_url])
        ch_db_props = ch_db.map { dir, url ->
            def props = dir.listFiles().find { it.name.endsWith('.properties') }
            if (!props) error "No .properties file found in rdp_db_dir: ${dir}"
            props
        }
    } else {
        def db_url = (params.its_region == 'ITS1') ? params.rdp_db_its1_url : params.rdp_db_its2_url
        DOWNLOAD_DB(
            Channel.value(db_url),
            Channel.value(params.its_region)
        )
        ch_db_props = DOWNLOAD_DB.out.db_dir.map { db_dir ->
            def props = db_dir.listFiles().find { it.name.endsWith('.properties') }
            if (!props) error "No .properties file found in downloaded db dir: ${db_dir}"
            props
        }
        // versions.yml written inside db_dir by DOWNLOAD_DB
    }

    // ── Per-sample trimming ───────────────────────────────────────────────
    TRIM_GALORE(ch_reads)
    ch_versions = ch_versions.mix(TRIM_GALORE.out.versions.first())

    // ── Joint DADA2 denoising ─────────────────────────────────────────────
    ch_r1_all = TRIM_GALORE.out.reads.map { meta, r1, r2 -> r1 }.collect()
    ch_r2_all = TRIM_GALORE.out.reads.map { meta, r1, r2 -> r2 }.collect()

    DADA2(ch_r1_all, ch_r2_all)
    ch_versions = ch_versions.mix(DADA2.out.versions)

    // ── ITS extraction ────────────────────────────────────────────────────
    ITSX(
        DADA2.out.asvs,
        params.its_region,
        params.its_min_len,
        params.its_max_len
    )
    ch_versions = ch_versions.mix(ITSX.out.versions)

    // ── Consolidate (ITSxed + not-ITSxed) ────────────────────────────────
    CONSOLIDATE(DADA2.out.asvs, ITSX.out.itsxed)
    ch_versions = ch_versions.mix(CONSOLIDATE.out.versions)

    // ── RDP classification ────────────────────────────────────────────────
    RDP_CLASSIFIER(
        CONSOLIDATE.out.fasta,
        ch_db_props,
        params.rdp_confidence
    )
    ch_versions = ch_versions.mix(RDP_CLASSIFIER.out.versions)

    // ── Filter count table to classified ASVs ─────────────────────────────
    FILTER_ASV_TABLE(
        DADA2.out.counts,
        CONSOLIDATE.out.fasta
    )
    ch_versions = ch_versions.mix(FILTER_ASV_TABLE.out.versions)

    // ── Merge ASVs with identical taxonomy ───────────────────────────────
    MERGE_DUPLICATES(
        FILTER_ASV_TABLE.out.counts,
        RDP_CLASSIFIER.out.taxonomy
    )
    ch_versions = ch_versions.mix(MERGE_DUPLICATES.out.versions)

    emit:
    asvs          = CONSOLIDATE.out.fasta
    counts        = FILTER_ASV_TABLE.out.counts
    taxonomy      = RDP_CLASSIFIER.out.taxonomy
    merged_counts = MERGE_DUPLICATES.out.counts
    trim_stats    = TRIM_GALORE.out.stats
    dada2_stats   = DADA2.out.stats
    itsx_summary  = ITSX.out.summary
    versions      = ch_versions
}
