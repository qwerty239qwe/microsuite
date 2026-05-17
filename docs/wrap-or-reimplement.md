# Wrap or reimplement policy

Version 0.1.0 reimplements small, stable methods where correctness is easy to
test:

- observed features
- Shannon diversity
- Simpson diversity
- Pielou evenness
- Bray-Curtis distance
- Jaccard distance

Heavy, platform-sensitive, or fast-moving methods are wrapped instead:

- BIOM import uses `biom-format`
- QIIME2 `.qza` import reads artifact layouts and delegates BIOM tables when needed
- ANCOM-BC is executed through external `Rscript`

Future releases can port selected methods after the CLI, data model, and tests
are stable.
