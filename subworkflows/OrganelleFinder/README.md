# OrganelleFinder pipeline

## Quickstart with docker

```bash
nextflow run -profile docker /path/to/OrganelleFinder.nf \
  -params-file '/path/to/params_file.yml'
```

## Usage

### Parameters

- General:
   * `genome_assembly`: Path to FNA genome assembly.
   * `reference_mitochondria`: Path to FNA reference mitochondria.
   * `reference_chloroplast`: Path to FNA reference chloroplast (if plant).
   * `reads_file`: Path to PacBio reads file (if present).
   * `input_type`: Type of organism, either 'animal' or 'plant'.

- Mitochondrial parameters:
   * `mit_blast_evalue`: E-value threshold for the mitochondria matches (default: 1e-6).
   * `mit_bitscore`: Bit score threshold for the mitochondria matches (default plant: 100, default animal: 100).
   * `mit_significant_gene_matches`: Threshold for number of unique matches for contig to be classified as mitochondrial (default plant: 5, default animal: 2).
   * `mit_suspicious_gene_matches`: Threshold for number of unique matches for contig to be classified as suspicious (default plant: 2, default animal: 1).
   * `mit_max_contig_length`: Max bp length of contig to be classified as mitochondrial (default plant: 2500000, default animal: 100000)
   * `mit_min_span_fraction`: Minimal match span fraction for contig to be classified as mitochondrial (default plant: 0.8, default animal: 0.8)

- Chloroplast parameters:
   * `chl_blast_evalue`: E-value threshold for the chloroplast matches (default: 1e-6).
   * `chl_bitscore`: Bit score threshold for the chloroplast matches (default: 150).
   * `chl_significant_gene_matches`: Threshold for number of unique matches for contig to be classified as chloroplast (default: 9).
   * `chl_suspicious_gene_matches`: Threshold for number of unique matches for contig to be classified as suspicious (default: 2).
   * `chl_max_contig_length`: Max bp length of contig to be classified as chloroplast (default: 400000)
   * `chl_min_span_fraction`: Minimal match span fraction for contig to be classified as chloroplast (default: 0.8)


Parameters to the workflow is provided by via a params-file that may look as follows:

`animal_params.yml`:
```
// General parameters
genome_assembly : ''
reference_mitochondria : ''
reads_file : ''
input_type : 'animal'

// Mitochondrial parameters
mit_blast_evalue : '1e-6'
mit_bitscore : 100
mit_significant_gene_matches : 2
mit_suspicious_gene_matches : 1
mit_max_contig_length : 100000
mit_min_span_fraction : 0.8

```

