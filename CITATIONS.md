# repgenr: Citations

## Pipeline

RepGenR builds representative-genome repositories and the associated FlexTaxD
taxonomy. If you use it, please cite this repository:
<https://github.com/druvus/repgenr>.

## Workflow managers

- [Nextflow](https://doi.org/10.1038/nbt.3820)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C.
  > Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017.

## Tools

The stages below shell out to external tools; cite the ones your run uses.

- **Genome selection / download**
  - [GTDB](https://doi.org/10.1093/nar/gkab776) (bacterial taxonomy / metadata)
  - [BV-BRC](https://doi.org/10.1093/nar/gkac1003) (viral metadata)
  - [NCBI Datasets](https://doi.org/10.1093/nar/gkad924) (genome download)

- **Dereplication**
  - [skDER / skani](https://doi.org/10.1186/s13059-024-03184-z)
  - [Galah](https://doi.org/10.1186/s40168-021-01213-8)
  - [sourmash](https://doi.org/10.21105/joss.00027) and the
    [branchwater plugin](https://doi.org/10.21105/joss.06830)
  - [dRep](https://doi.org/10.1038/ismej.2017.126)

- **Alignment and phylogeny**
  - [progressiveMauve](https://doi.org/10.1371/journal.pone.0011147)
  - [SibeliaZ](https://doi.org/10.1038/s41467-020-19777-8)
  - [Cactus / Minigraph-Cactus](https://doi.org/10.1038/s41587-023-01793-w)
  - [Mashtree](https://doi.org/10.21105/joss.01762)
  - [RAxML-NG](https://doi.org/10.1093/bioinformatics/btz305)
  - [IQ-TREE 2](https://doi.org/10.1093/molbev/msaa015)
  - [FastTree 2](https://doi.org/10.1371/journal.pone.0009490)

- **SNP typing**
  - [minimap2](https://doi.org/10.1093/bioinformatics/bty191),
    [SAMtools / BCFtools](https://doi.org/10.1093/gigascience/giab008)
  - [Snippy](https://github.com/tseemann/snippy)
  - [ParSNP / Harvest](https://doi.org/10.1186/s13059-014-0524-x)
  - [Gubbins](https://doi.org/10.1093/nar/gku1196) (recombination masking)

- **Taxonomy**
  - [FlexTaxD](https://doi.org/10.1093/bioinformatics/btab372)

## Software packaging

- [BioContainers](https://doi.org/10.1093/bioinformatics/btx192)
- [Seqera Wave](https://seqera.io/wave/)
- [Anaconda / conda-forge / Bioconda](https://doi.org/10.1038/s41592-018-0046-7)
