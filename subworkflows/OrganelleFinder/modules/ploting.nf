process PLOTING {

    conda "${task.ext.enable_conda  ? 'bioconda::tool=bioinfokit:2.0.8' : '' }"
    container "${workflow.containerEngine == 'singularity' &&
                  !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/bioinfokit:2.0.8--pyh5e36f6f_0' :
              'quay.io/biocontainers/bioinfokit:2.0.8--pyh5e36f6f_0' }"

    input:
    path alignment

    output:
    path "*.png*"

    script: 
    """
    #!/usr/bin/env python
    import sys
    import pandas as pd
    import matplotlib.pyplot as plt

    files = "${alignment}".split(" ")
    for file in files:
        image_name = file.replace(".tsv", ".png")
        csvFile = open(file, 'r')
        df = pd.read_csv(csvFile, sep='\t', names=["Accession", "Position", "Depth"])
        csvFile.close()
        fig, ax = plt.subplots()
        fig.set_size_inches(50, 5)
        df.plot(kind='scatter',x='Position',y='Depth',ax=ax)
        plt.ylabel('Coverage')
        plt.xlabel('Position')
        plt.title(image_name + ' coverage plot')
        plt.savefig(image_name)
    """
}
