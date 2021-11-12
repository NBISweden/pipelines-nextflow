#!/bin/bash -ue
blastx -query GCF_000002985.6_WBcel235_genomic.fna -db Ref_mitochondria_animal.fna Ref_mitochondria_animal.fna.pdb Ref_mitochondria_animal.fna.phr Ref_mitochondria_animal.fna.pin Ref_mitochondria_animal.fna.pot Ref_mitochondria_animal.fna.psq Ref_mitochondria_animal.fna.ptf Ref_mitochondria_animal.fna.pto  -out GCF_000002985.6_WBcel235_genomic_blast.tsv
