process MS_FUNCTIONAL {
    tag 'ms_functional'
    label 'microsuite_picrust2'
    cpus { params.threads as int }
    publishDir "${params.outdir}/functional", mode: 'copy'

    input:
    path otu_table
    path rep_seqs

    output:
    path 'picrust2', emit: functional_dir

    script:
    def shell_quote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }
    def enabled = { value ->
        value instanceof Boolean
            ? value
            : value != null && value.toString().trim().toLowerCase() in ['1', 'true', 'yes', 'on']
    }
    def argv = [
        '--backend', 'picrust2',
        '--table', otu_table,
        '--rep-seqs', rep_seqs,
        '--output-dir', 'picrust2',
        '--threads', task.cpus,
        '--picrust2-database', params.picrust2_database ?: 'SC',
    ]
    if (params.picrust2_max_nsti != null && params.picrust2_max_nsti.toString().trim()) {
        argv += ['--picrust2-max-nsti', params.picrust2_max_nsti]
    }
    if (enabled(params.picrust2_no_pathways)) {
        argv += ['--picrust2-no-pathways']
    }
    if (enabled(params.picrust2_coverage)) {
        argv += ['--picrust2-coverage']
    }
    if (enabled(params.picrust2_no_regroup)) {
        argv += ['--picrust2-no-regroup']
    }
    """
    # SC and oldIMG deliberately retain the historical two-input process call.
    microsuite functional_profile ${argv.collect(shell_quote).join(' ')}
    """

    stub:
    """
    mkdir -p picrust2
    printf 'stub picrust2 output\\n' > picrust2/pathways.tsv
    """
}

// SC/oldIMG mapping overrides use this process. It stages path-valued maps in
// the task directory while allowing reaction_func to remain an upstream
// symbolic value such as EC or KO.
process MS_FUNCTIONAL_OVERRIDES {
    tag 'ms_functional_overrides'
    label 'microsuite_picrust2'
    cpus { params.threads as int }
    publishDir "${params.outdir}/functional", mode: 'copy'

    input:
    path otu_table
    path rep_seqs
    path override_assets, stageAs: 'picrust2_override_asset??/*'
    val override_plan

    output:
    path 'picrust2', emit: functional_dir

    script:
    def shell_quote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }
    def enabled = { value ->
        value instanceof Boolean
            ? value
            : value != null && value.toString().trim().toLowerCase() in ['1', 'true', 'yes', 'on']
    }
    def staged_at = { index ->
        index == null ? null : override_assets[index as int]
    }
    def argv = [
        '--backend', 'picrust2',
        '--table', otu_table,
        '--rep-seqs', rep_seqs,
        '--output-dir', 'picrust2',
        '--threads', task.cpus,
        '--picrust2-database', params.picrust2_database ?: 'SC',
    ]
    if (override_plan.pathway_map != null) {
        argv += ['--picrust2-pathway-map', staged_at(override_plan.pathway_map)]
    }
    if (override_plan.reaction_func_path != null) {
        argv += ['--picrust2-reaction-func', staged_at(override_plan.reaction_func_path)]
    } else if (override_plan.reaction_func_value != null && override_plan.reaction_func_value.toString().trim()) {
        argv += ['--picrust2-reaction-func', override_plan.reaction_func_value]
    }
    if (override_plan.regroup_map != null) {
        argv += ['--picrust2-regroup-map', staged_at(override_plan.regroup_map)]
    }
    if (params.picrust2_max_nsti != null && params.picrust2_max_nsti.toString().trim()) {
        argv += ['--picrust2-max-nsti', params.picrust2_max_nsti]
    }
    if (enabled(params.picrust2_no_pathways)) {
        argv += ['--picrust2-no-pathways']
    }
    if (enabled(params.picrust2_coverage)) {
        argv += ['--picrust2-coverage']
    }
    if (enabled(params.picrust2_no_regroup)) {
        argv += ['--picrust2-no-regroup']
    }
    """
    # Mapping overrides are rewritten to staged paths before shell quoting.
    microsuite functional_profile ${argv.collect(shell_quote).join(' ')}
    """

    stub:
    """
    mkdir -p picrust2
    printf 'stub picrust2 output\\n' > picrust2/pathways.tsv
    """
}

// A symbolic reaction_func does not need an input file. Keep this process
// separate so an empty path collection is never required for a valid run.
process MS_FUNCTIONAL_SYMBOLIC_OVERRIDE {
    tag 'ms_functional_symbolic_override'
    label 'microsuite_picrust2'
    cpus { params.threads as int }
    publishDir "${params.outdir}/functional", mode: 'copy'

    input:
    path otu_table
    path rep_seqs

    output:
    path 'picrust2', emit: functional_dir

    script:
    def shell_quote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }
    def enabled = { value ->
        value instanceof Boolean
            ? value
            : value != null && value.toString().trim().toLowerCase() in ['1', 'true', 'yes', 'on']
    }
    def argv = [
        '--backend', 'picrust2',
        '--table', otu_table,
        '--rep-seqs', rep_seqs,
        '--output-dir', 'picrust2',
        '--threads', task.cpus,
        '--picrust2-database', params.picrust2_database ?: 'SC',
        '--picrust2-reaction-func', params.picrust2_reaction_func,
    ]
    if (params.picrust2_max_nsti != null && params.picrust2_max_nsti.toString().trim()) {
        argv += ['--picrust2-max-nsti', params.picrust2_max_nsti]
    }
    if (enabled(params.picrust2_no_pathways)) {
        argv += ['--picrust2-no-pathways']
    }
    if (enabled(params.picrust2_coverage)) {
        argv += ['--picrust2-coverage']
    }
    if (enabled(params.picrust2_no_regroup)) {
        argv += ['--picrust2-no-regroup']
    }
    """
    # reaction_func is restricted to a symbolic value before this process is
    # selected, so no host path can be interpolated here.
    microsuite functional_profile ${argv.collect(shell_quote).join(' ')}
    """

    stub:
    """
    mkdir -p picrust2
    printf 'stub picrust2 output\\n' > picrust2/pathways.tsv
    """
}

// Custom references are a separate process because every custom asset must be
// a real Nextflow path input. The value plan contains only roles and indexes;
// the command below resolves those indexes against the staged path collection.
// This avoids interpolating host paths into a Docker/Singularity task.
process MS_FUNCTIONAL_CUSTOM {
    tag 'ms_functional_custom'
    label 'microsuite_picrust2'
    cpus { params.threads as int }
    publishDir "${params.outdir}/functional", mode: 'copy'

    input:
    path otu_table
    path rep_seqs
    // Use one indexed subdirectory per asset so duplicate filenames from two
    // custom domains cannot collide during staging.
    path custom_assets, stageAs: 'picrust2_custom_asset??/*'
    val custom_plan

    output:
    path 'picrust2', emit: functional_dir

    script:
    def shell_quote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }
    def enabled = { value ->
        value instanceof Boolean
            ? value
            : value != null && value.toString().trim().toLowerCase() in ['1', 'true', 'yes', 'on']
    }
    def staged_at = { index ->
        index == null ? null : custom_assets[index as int]
    }
    def argv = [
        '--backend', 'picrust2',
        '--table', otu_table,
        '--rep-seqs', rep_seqs,
        '--output-dir', 'picrust2',
        '--threads', task.cpus,
        '--picrust2-database', 'custom',
        '--picrust2-ref-dir1', staged_at(custom_plan.ref_dir1),
    ]
    if (custom_plan.ref_dir2 != null) {
        argv += ['--picrust2-ref-dir2', staged_at(custom_plan.ref_dir2)]
    }
    (custom_plan.trait_tables_ref1 ?: []).each { index ->
        argv += ['--picrust2-custom-trait-tables-ref1', staged_at(index)]
    }
    (custom_plan.trait_tables_ref2 ?: []).each { index ->
        argv += ['--picrust2-custom-trait-tables-ref2', staged_at(index)]
    }
    argv += ['--picrust2-marker-gene-table-ref1', staged_at(custom_plan.marker_gene_table_ref1)]
    if (custom_plan.marker_gene_table_ref2 != null) {
        argv += ['--picrust2-marker-gene-table-ref2', staged_at(custom_plan.marker_gene_table_ref2)]
    }
    if (custom_plan.pathway_map != null) {
        argv += ['--picrust2-pathway-map', staged_at(custom_plan.pathway_map)]
    }
    if (custom_plan.reaction_func_path != null) {
        argv += ['--picrust2-reaction-func', staged_at(custom_plan.reaction_func_path)]
    } else if (custom_plan.reaction_func_value != null && custom_plan.reaction_func_value.toString().trim()) {
        argv += ['--picrust2-reaction-func', custom_plan.reaction_func_value]
    }
    if (custom_plan.regroup_map != null) {
        argv += ['--picrust2-regroup-map', staged_at(custom_plan.regroup_map)]
    }
    if (params.picrust2_max_nsti != null && params.picrust2_max_nsti.toString().trim()) {
        argv += ['--picrust2-max-nsti', params.picrust2_max_nsti]
    }
    if (enabled(params.picrust2_no_pathways)) {
        argv += ['--picrust2-no-pathways']
    }
    if (enabled(params.picrust2_coverage)) {
        argv += ['--picrust2-coverage']
    }
    if (enabled(params.picrust2_no_regroup)) {
        argv += ['--picrust2-no-regroup']
    }
    """
    # Every custom path in custom_assets is staged by Nextflow before this
    # command runs; staged_at() rewrites the option values to those paths.
    microsuite functional_profile ${argv.collect(shell_quote).join(' ')}
    """

    stub:
    """
    mkdir -p picrust2
    printf 'stub picrust2 output\\n' > picrust2/pathways.tsv
    """
}
