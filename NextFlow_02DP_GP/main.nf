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
    conda '~/.conda/envs/DPGP_env'

    input:
        file gem from input_gems

    script:
    """
    iteration=1000

    python \${DP_GP_PATH}/DP_GP_cluster.py\
            -i \${gem} \
            -o \${gem.baseName} \
            -n \${iteration} \
            --fast \
            --true_times
    """
}