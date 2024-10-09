# Biothings API

## Overview

[BioThings.io](https://biothings.io) is a collection of APIs that provide comprehensive data annotation services for various types of biomedical information. It is designed to facilitate data access and integration for researchers working in genomics, drug discovery, and related fields. BioThings APIs are optimized for high-speed queries and data annotation, making it an essential tool for accessing gene, variant, chemical, and disease-related information.

### Key BioThings APIs

BioThings offers four main services, each specializing in a specific type of annotation:

1. **[MyGene.info](https://mygene.info/)**: Provides gene-centric annotation, including gene descriptions, biological functions, interactions, and pathway information.
2. **[MyVariant.info](https://myvariant.info/)**: Focuses on variant annotation, offering detailed insights into genetic variations, their impacts, and clinical significance.
3. **[MyChem.info](https://mychem.info/)**: A chemical and drug annotation service that provides data on chemical structures, drug interactions, and pharmacological information.
4. **[MyDisease.info](https://mydisease.info/)**: Specializes in disease annotation, linking diseases to associated genes, phenotypes, clinical trials, and more.

Each of these services is designed to be scalable, fast, and accessible via a simple RESTful API interface.

## Services Overview

### 1. MyGene.info

- **Description**: MyGene.info offers a powerful gene annotation service, providing researchers access to comprehensive gene-centric data. It aggregates information from various sources such as NCBI, Ensembl, UniProt, and more.
- **Key Features**:
  - Gene descriptions and summaries
  - Biological functions and interactions
  - Pathway information
  - Support for human and model organism genes

### 2. MyVariant.info

- **Description**: MyVariant.info is a variant-centric service that delivers high-quality annotations for genetic variations. It aggregates variant data from multiple sources, including dbSNP, ClinVar, and ExAC.
- **Key Features**:
  - Population frequencies
  - Clinical significance and pathogenicity
  - Functional predictions and annotations
  - Support for SNPs, indels, and structural variants

### 3. MyChem.info

- **Description**: MyChem.info provides an efficient way to retrieve chemical and drug annotations. It integrates data from ChEMBL, DrugBank, PubChem, and other chemical databases.
- **Key Features**:
  - Chemical structure and formula
  - Drug interactions and adverse effects
  - Pharmacological properties and indications
  - Toxicology and safety information

### 4. MyDisease.info

- **Description**: MyDisease.info is a disease-centric annotation service designed to provide detailed information about human diseases, including gene-disease associations, clinical trials, and ontology information.
- **Key Features**:
  - Disease descriptions and ontology terms
  - Gene and variant associations
  - Links to clinical trials and disease models
  - Disease-gene interaction networks

## How to Access

Each of the BioThings APIs provides a well-documented RESTful interface, allowing users to send simple HTTP queries to retrieve data in JSON format. The APIs support a wide range of query options, including gene/variant/disease names, chemical identifiers, and specific properties.

Biothing's API (https://biothings.io/)

MyGene.info (https://mygene.info/)

MyVariant.info (https://myvariant.info/)

MyChem.info (https://mychem.info/)

MyDisease.info (https://mydisease.info/)



## Querying Genes

You can query genes from BioThings' MyGene.info API using a variety of identifiers and attributes. Here are the main types of queries you can run:

### 1. **Retrieve Gene Information by Specific ID**

   - You can look up detailed gene information using various gene identifiers:
     - **Entrez Gene ID**: `https://mygene.info/v3/gene/1017`
     - **Ensembl Gene ID**: Similar format to Entrez ID queries.
     - **UniProt ID**: Example query: `https://mygene.info/v3/query?q=uniprot:P24941`

### 2. **Search by Gene Symbol**

   - Retrieve genes by their symbol (e.g., `CDK2`):  
     - `https://mygene.info/v3/query?q=symbol:cdk2`

### 3. **Search by Gene Summary Text**

   - Search for genes based on their descriptions or related summaries. For example, to get genes related to "Insulin":
     - `https://mygene.info/v3/query?q=summary:insulin`

### 4. **Retrieve Data Source Metadata**

   - Retrieve metadata about the sources of gene data:
     - `https://mygene.info/metadata`

### 5. **Get a List of All Available Fields**

   - Get a list of all fields available in the MyGene.info database:
     - `https://mygene.info/metadata/fields`

These queries provide flexible options for retrieving gene-centric data using different identifiers, symbols, or biological descriptions.



## Querying Variants

You can query variants from BioThings' MyVariant.info API using various types of identifiers, genomic coordinates, and annotation criteria. Below is a summary of the main types of queries:

### 1. **Retrieve Variant by Specific HGVS ID**

   - Look up variant information using a specific HGVS (Human Genome Variation Society) ID:
     - `https://myvariant.info/v1/variant/chr7:g.55241707G>T`

### 2. **Retrieve Variants Based on Specific Annotations**

   - Query for variants that have specific annotations, such as those annotated in ClinVar:
     - `https://myvariant.info/v1/query/?q=_exists_:clinvar`

### 3. **Find Variants for a Specific Gene**

   - Retrieve all non-synonymous variants for a given gene using its symbol, e.g., `BTK`:
     - `https://myvariant.info/v1/query/?q=dbnsfp.genename:BTK`

### 4. **Search Variants by Genomic Range**

   - Search for variants within a specified genomic range:
     - `https://myvariant.info/v1/query/?q=chr1:69000-70000`

### 5. **Complex Queries Combining Multiple Annotations**

   - Create advanced queries that combine multiple annotations (e.g., variants with `wellderly` annotation and `possibly damaging` CADD PolyPhen predictions):
     - `https://myvariant.info/v1/query/?q=_exists_:wellderly AND cadd.polyphen.cat:possibly_damaging&fields=wellderly,cadd.polyphen`

### 6. **Retrieve Data Source Metadata**

   - Retrieve metadata about the sources of variant data:
     - `https://mygene.info/metadata/`

### 7. **Get a List of All Available Fields**

   - Get a complete list of all fields available in MyVariant.info:
     - `https://mygene.info/metadata/fields`

These queries enable a comprehensive exploration of variant data using various identifiers, genomic ranges, annotations, and gene symbols.



## Example Gene querying results

Here's a categorized table summarizing the useful information from the given JSON file:

| **Category**                  | **Field Name**     | **Description/Content**                                      |
| ----------------------------- | ------------------ | ------------------------------------------------------------ |
| **Identifiers**               | `AllianceGenome`   | 1771                                                         |
|                               | `HGNC`             | 1771                                                         |
|                               | `MIM`              | 116953                                                       |
|                               | `_id`              | 1017                                                         |
|                               | `entrezgene`       | 1017                                                         |
|                               | `ensembl`          | Gene: ENSG00000123374, Protein: multiple ENSP, Transcript: multiple ENST |
|                               | `symbol`           | CDK2 (Cyclin Dependent Kinase 2)                             |
|                               | `other_names`      | cdc2-related protein kinase, cell division protein kinase 2, p33 protein kinase |
| **Accession Numbers**         | `accession`        | Genomic, Protein, RNA, Translation                           |
|                               | `protein`          | Various accession numbers: e.g., AAA35667.1, NP_001277159.1  |
|                               | `rna`              | Various accession numbers: e.g., NM_001290230.2, XM_011537732.2 |
| **Aliases**                   | `alias`            | CDKN2, p33(CDK2)                                             |
| **Chemical Data**             | `chembl`           | Multiple entries linking to ChEMBL data, e.g., CHEMBL4106153 |
| **Genomic Position**          | `genomic_pos`      | chr 12: 55966781-55972789, Strand: 1                         |
|                               | `genomic_pos_hg19` | chr 12: 56360553-56366568, Strand: 1                         |
| **Ensembl Information**       | `ensembl`          | Gene: ENSG00000123374, multiple protein and transcript identifiers |
|                               | `type_of_gene`     | protein_coding                                               |
| **Functional Annotations**    | `go`               | Gene Ontology (GO) annotations for BP, MF, and CC categories |
|                               | `BP`               | Biological processes like G1/S transition, protein phosphorylation, signal transduction |
|                               | `CC`               | Cellular components: nucleus, centrosome, cyclin-CDK complex |
|                               | `MF`               | Molecular functions: ATP binding, cyclin binding, protein kinase activity |
| **Orthologs**                 | `orthologs`        | Homologous genes across multiple species                     |
| **Pathway Information**       | `pathway`          | Various pathway identifiers: BioCarta, KEGG, NetPath, PID, Reactome, WikiPathways |
| **Protein Domains**           | `interpro`         | Protein domain information, e.g., Serine/threonine-protein kinase, ATP binding site |
| **Protein Structure**         | `pdb`              | Multiple PDB structure identifiers for CDK2                  |
| **Pharmaceutical Data**       | `pharmgkb`         | PA101                                                        |
| **Proteomics**                | `refseq`           | Reference sequence identifiers for genomic, protein, and RNA |
|                               | `translation`      | Protein to RNA translation mappings                          |
| **Gene Information**          | `summary`          | Brief gene description                                       |
|                               | `map_location`     | 12q13.2                                                      |
| **Protein Interactions**      | `ipi`              | IPI identifiers                                              |
| **Transcript Information**    | `exons`            | Genomic exon positions for various transcripts               |
|                               | `exons_hg19`       | Exon positions in HG19 genome version                        |
|                               | `transcript`       | Associated transcripts, e.g., ENST00000266970.4              |
| **Gene Function Information** | `generif`          | PubMed-based gene references for functional studies          |
|                               | `pubmed`           | PubMed IDs                                                   |
|                               | `text`             | Description of study findings                                |
| **Evolutionary Data**         | `homologene`       | Homologene gene entries                                      |
|                               | `id`               | 74409                                                        |
|                               | `genes`            | Multiple homologous genes                                    |
| **Protein Family**            | `pantherdb`        | PANTHER family information                                   |
| **Experimental Data**         | `reagent`          | Experimental reagents (siRNA, shRNA, etc.) used for the gene |
| **Wikipedia Information**     | `wikipedia`        | URL: `Cyclin-dependent kinase 2`                             |
| **Protein Interactions**      | `reporter`         | Gene expression array probe identifiers                      |
| **UniProt Data**              | `uniprot`          | Swiss-Prot: P24941, TrEMBL: multiple entries                 |
|                               | `unii`             | IEQ8K30M47                                                   |



## Example Variant querying results

Here is a table summarizing the key categories and fields from the given JSON file:

| **Category**                  | **Field**                              | **Value/Description**                                        |
| ----------------------------- | -------------------------------------- | ------------------------------------------------------------ |
| **General Information**       | `_id`                                  | `chr7:g.55241707G>T`                                         |
|                               | `_version`                             | `2`                                                          |
| **Variant Type**              | `type`                                 | `SNV` (Single Nucleotide Variant)                            |
|                               | `ref`                                  | `G`                                                          |
|                               | `alt`                                  | `T`                                                          |
|                               | `chrom`                                | `7`                                                          |
|                               | `pos`                                  | `55241707`                                                   |
|                               | `rsid`                                 | `rs28929495`                                                 |
|                               | `observed`                             | `true`                                                       |
| **CADD Scores**               | `cadd.phred`                           | `26.2`                                                       |
|                               | `cadd.annotype`                        | `CodingTranscript`                                           |
|                               | `cadd.bstatistic`                      | `822`                                                        |
|                               | `cadd.consdetail`                      | `missense`                                                   |
|                               | `cadd.consequence`                     | `NON_SYNONYMOUS`                                             |
|                               | `cadd.gc`                              | `0.57`                                                       |
| **Gene Information**          | `cadd.gene.genename`                   | `EGFR`                                                       |
|                               | `cadd.gene.ccds_id`                    | `CCDS5514.1`                                                 |
|                               | `cadd.gene.feature_id`                 | `ENST00000275493`                                            |
|                               | `cadd.gene.gene_id`                    | `ENSG00000146648`                                            |
|                               | `cadd.gene.cds.cds_pos`                | `2155`                                                       |
|                               | `cadd.gene.prot.protpos`               | `719`                                                        |
| **Protein Change**            | `clinvar.hgvs.protein`                 | Various entries, e.g., `NP_001333827.1:p.Gly719Cys`          |
|                               | `cadd.oaa`                             | `G` (Original Amino Acid)                                    |
|                               | `cadd.naa`                             | `C` (New Amino Acid)                                         |
|                               | `cgi.protein_change`                   | `EGFR:G719C`                                                 |
| **Clinical Information**      | `clinvar.rcv[].clinical_significance`  | Various entries: `drug response`, `Pathogenic`, `Likely pathogenic` |
|                               | `clinvar.rcv[].conditions.name`        | Various entries: `Nonsmall cell lung cancer, response to tyrosine kinase inhibitor in, somatic` |
|                               | `clinvar.rcv[].origin`                 | `somatic`                                                    |
| **Drug Response**             | `cgi[].association`                    | `Responsive`                                                 |
|                               | `cgi[].drug`                           | Various entries: `HSP90 inhibitors`, `EGFR TK inhibitors`, `Afatinib`, `Gefitinib` |
|                               | `cgi[].primary_tumor_type`             | Various entries: `Lung`, `Non-small cell lung`               |
|                               | `cgi[].evidence_level`                 | Various entries: `FDA guidelines`, `Early trials`, `NCCN guidelines` |
| **Functional Annotations**    | `dbnsfp.aa.alt`                        | `C`                                                          |
|                               | `dbnsfp.appris`                        | `principal1`                                                 |
|                               | `dbnsfp.eigen.raw_coding`              | `0.953480363944939`                                          |
| **Pathogenicity Predictions** | `polyphen.cat`                         | `probably_damaging`                                          |
|                               | `polyphen.val`                         | `0.998`                                                      |
|                               | `sift.cat`                             | `deleterious`                                                |
|                               | `sift.val`                             | `0`                                                          |
| **Conservation Scores**       | `gerp.rs`                              | `387.1`                                                      |
|                               | `phylop.mammalian`                     | `2.742`                                                      |
|                               | `phast_cons.mammalian`                 | `1.0`                                                        |
|                               | `fitcons`                              | `0.604399`                                                   |
| **Regulatory Annotations**    | `cadd.encode.h3k27ac`                  | `2.48`                                                       |
|                               | `cadd.encode.exp`                      | `341.62`                                                     |
|                               | `cadd.encode.nucleo`                   | `1.8`                                                        |
|                               | `cadd.mapability.20bp`                 | `1`                                                          |
| **Pharmacogenomics**          | `clinvar.rcv[].conditions.identifiers` | Various entries, e.g., `medgen: C4016032`, `mondo: MONDO:0005233` |
| **Publications**              | `dbsnp.citations`                      | Various PubMed IDs: `15118073, 15329413`                     |
|                               | `docm.pubmed_id`                       | `15118073, 19922469`                                         |
|                               | `mutdb.cosmic_id`                      | `20881, 6253`                                                |
|                               | `clinvar.variant_id`                   | `16611`                                                      |



## Documentation and Community

For detailed documentation and usage examples, visit [BioThings.io Documentation](https://biothings.io/explorer). The community is active and provides support through forums, tutorials, and GitHub repositories.

## Conclusion

BioThings.io offers a comprehensive suite of annotation services for researchers and bioinformaticians, enabling easy access to high-quality data for genes, variants, chemicals, and diseases. Its APIs are built to be fast, scalable, and easy to integrate, making it a valuable resource for a wide range of biomedical research applications.