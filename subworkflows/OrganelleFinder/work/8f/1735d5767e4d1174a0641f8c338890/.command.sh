#!/bin/bash -ue
blastx -query GCF_000002985.6_WBcel235_genomic.fna -db Ref_mitochondria_animal.fna -out GCF_000002985.6_WBcel235_genomic_blast.tsv
