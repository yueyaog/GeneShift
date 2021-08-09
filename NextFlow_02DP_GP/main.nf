println """\
====================================
G E N E S H I F T   W O R K F L O W
====================================

Workflow Information:
---------------------
  Launch Directory:   ${workflow.launchDir}
  Work Directory:     ${workflow.workDir}
  Config Files:       ${workflow.configFiles}
  Profiles:           ${workflow.profile}


"""


/**
* Load GEM file for each DP_GP_cluster input
*/

input_gems = Channel.fromPath("${params.input.dir}/${params.input.gem_files}")


/**
 * The DP_GP process performs a single run of DP_GP_cluster
 * for each input GEM.
 */

process DP_GP {
    publishDir "${params.output.dir}"

    input:
    file gem from input_gems

    script:
    """
    module load anaconda3/5.1.0-gcc 
    source activate ${params.conda_env}
    python ${DPGP_PATH}/DP_GP_cluster.py\
            -i ${gem} \
            -o ${gem.baseName} \
            -n ${params.iteration} \
            --fast \
            --true_times
    """
}