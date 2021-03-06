workflow SctkMastDge {
    
    # An integer matrix of counts
    File raw_counts

    # The samples/cells in one group (required)
    Array[String] expSamples

    # The samples/cells in the other group.
    # Not required. If missing, then runs a one-vs-all
    Array[String]? baseSamples

    # a name assigned to the "experimental" group. Required
    String expGroupName

    # If we are performing an explicit group comparison, 
    # we require this.
    String? baseGroupName

    Boolean isDirectComparison = defined(baseSamples)

    call runMast {
        input:
            raw_counts = raw_counts,
            expSamples = expSamples,
            baseSamples = baseSamples,
            expGroupName = expGroupName,
            baseGroupName = baseGroupName,
            isDirectComparison = isDirectComparison
    }

    output {
        File mast_output = runMast.fout
    }
}

task runMast {
    File raw_counts
    Array[String] expSamples
    Array[String]? baseSamples
    String expGroupName
    String? baseGroupName
    Boolean isDirectComparison

    String output_name_prefix = "sctk_mast_results"
    Int disk_size = 20

    command {
        if [[ "${isDirectComparison}" == "true" ]]
        then
          Rscript /opt/software/mast_dge.R \
            -f ${raw_counts} \
            -o ${output_name_prefix} \
            -a ${sep="," expSamples} \
            -b ${sep="," baseSamples} \
            --experimental_group_name ${expGroupName} \
            ${"--base_group_name " + baseGroupName}
        else
          Rscript /opt/software/mast_dge.R \
            -f ${raw_counts} \
            -o ${output_name_prefix} \
            -a ${sep="," expSamples} \
            --experimental_group_name ${expGroupName} \
            ${"--base_group_name " + baseGroupName}

        fi
    }

    output {
        File fout = glob("${output_name_prefix}*")[0]
    }

    runtime {
        docker: "ghcr.io/web-mev/mev-sctk-mast"
        cpu: 2
        memory: "16 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
