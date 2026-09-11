# RepGenR documentation

| Page | Read it when |
|------|--------------|
| [usage.md](usage.md) | You want to run the pipeline: the CLI stages, starting from local genomes, viruses, resume and `--force`, representative selection, SNP typing and masking, tree2tax collapsing, the Nextflow layer (parameters, profiles, scaling, per-process configuration), running tools in containers, and troubleshooting. |
| [cli-reference.md](cli-reference.md) | You need the exact options of a command. Generated from the command tree. |
| [output.md](output.md) | You want to know what a stage wrote and what each file means. |
| [developing.md](developing.md) | You are changing the code: architecture, data contracts, adding a tool adapter. |
| [verification.md](verification.md) | You want to know which adapters have been run against their real tools, on what data, and how long the measured runs took. |
| [audit/](audit/README.md) | Records behind the current design: the CLI matrix and the 2026 scaling and bias audit. |

The changelog is at the repository root ([CHANGELOG.md](../CHANGELOG.md)).
