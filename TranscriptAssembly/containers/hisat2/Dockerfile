FROM continuumio/miniconda3:4.8.2

SHELL ["/bin/bash", "-c"]

LABEL description="Hisat2 and Samtools container " \
      author="Mahesh Binzer-Panchal" 

RUN apt-get update --fix-missing && \
    apt-get install -y procps

RUN conda update -n base conda && \
    conda install -c conda-forge -c bioconda hisat2=2.1.0 samtools=1.9 && \
    conda clean --all -f -y

CMD [ "hisat2" ]
