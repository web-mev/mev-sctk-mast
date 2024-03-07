// default params. Can be overridden via the params JSON or commandline
params.baseSamples = []
params.baseGroupName = "ctrl"

process run_mast {

    tag "Run SCTK MAST tool"
    publishDir "${params.output_dir}/SctkMastDge.mast_output", mode:"copy", pattern:"${output_name_prefix}*"
    container "ghcr.io/web-mev/mev-sctk-mast"
    cpus 2
    memory '16 GB'

    input:
        path raw_counts

    output:
        path "${output_name_prefix}*"

    script:
      def base_samples_size = params.baseSamples.size()
      def flag = ""
      if ( base_samples_size > 0 ){
        flag = "-b ${params.baseSamples.join(',')}"
      }
      output_name_prefix = "sctk_mast_results"
        """
        Rscript /opt/software/mast_dge.R \
          -f ${raw_counts} \
          -o ${output_name_prefix} \
          -a ${params.expSamples.join(',')} \
          ${flag} \
          --experimental_group_name ${params.expGroupName} \
          ${"--base_group_name " + params.baseGroupName}
        """
}

workflow {
    run_mast(params.raw_counts)
}