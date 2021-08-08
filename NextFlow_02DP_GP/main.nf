println """\
=================================
 K I N C - N F   P I P E L I N E
=================================

Workflow Information:
---------------------
  Launch Directory:   ${workflow.launchDir}
  Work Directory:     ${workflow.workDir}
  Config Files:       ${workflow.configFiles}
  Profiles:           ${workflow.profile}

Execution Parameters:
---------------------
input:
    GEM:    ${params.input.GEMs}
"""


/**
* Load GEM file for each DP_GP_cluster input
*/
GEM_FILES = Channel.fromFilePairs("${params.input.dir}/${params.input.gem_files}", size: 1, flat: true)

/**
 * The DP_GP process performs a single run of DP_GP_cluster
 * for each input GEM.
 */

process DP_GP {
    publishDir "${params.output.dir}"

    input:
        each(c) from CONDITIONS
        set val(geometry), file(gmy_file), file(xml_file) from DATASETS
        each(trial) from Channel.from( 0 .. params.input.trials-1 )

    script:
    """
    input_GEM=\${BASEDIR}/02-DP_GP/Inputs/KValDPGP_input/clusterX.txt
    output_file=\${BASEDIR}/02-DP_GP/KValDPGP_output/clusterX
    iteration=1000

    python \${DP_GP_PATH}/DP_GP_cluster.py\
            -i \${input_GEM} \
            -o \${output_file} \
            -n \${iteration} \
            --fast \
            --true_times
    """