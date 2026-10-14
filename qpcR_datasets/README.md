# Bundled qPCR datasets

The files in this directory are benchmark qPCR datasets redistributed from the R package
[**qpcR**](https://cran.r-project.org/package=qpcR) by Andrej-Nikolai Spiess and Christian Ritz:

> Ritz, C. & Spiess, A.-N. (2008). qpcR: an R package for sigmoidal model selection in quantitative
> real-time polymerase chain reaction analysis. *Bioinformatics* 24, 1549–1551.

## License

These data files are distributed under the license of qpcR: **GPL (≥ 2)**, i.e. GPL-2 or GPL-3.
They are **not** covered by the MIT license that applies to the source code of qPCRModels.jl
(see [`../LICENSE`](../LICENSE)). The full license texts are available at
<https://www.gnu.org/licenses/old-licenses/gpl-2.0.html> and <https://www.gnu.org/licenses/gpl-3.0.html>.

## Sources

Experimental setup, instrument and assay details for every dataset are available in Julia through
`DatasetDetails["<name>"]`. The original publication of each dataset is given on the corresponding
qpcR help page (e.g. `?batsch` in R). The datasets taken from published studies include:

| Datasets | Original publication |
|---|---|
| `batsch1`–`batsch5` | Batsch, A. et al. (2008). Simultaneous fitting of real-time PCR data with efficiency of amplification modeled as Gaussian function of target fluorescence. *BMC Bioinformatics* 9, 95. |
| `boggy` | Boggy, G. J. & Woolf, P. J. (2010). A mechanistic model of PCR for accurate quantification of quantitative PCR data. *PLoS ONE* 5, e12355. |
| `guescini1`, `guescini2` | Guescini, M. et al. (2008). A new real-time PCR method to overcome significant quantitative inaccuracy due to slight amplification inhibition. *BMC Bioinformatics* 9, 326. |
| `karlen1`–`karlen3` | Karlen, Y. et al. (2007). Statistical significance of quantitative PCR. *BMC Bioinformatics* 8, 131. |
| `lievens1`–`lievens3` | Lievens, A. et al. (2012). Enhanced analysis of real-time PCR data by using a variable efficiency model: FPK-PCR. *Nucleic Acids Research* 40, e10. |
| `rutledge` | Rutledge, R. G. (2004). Sigmoidal curve-fitting redefines quantitative real-time PCR with the prospective of developing automated high-throughput applications. *Nucleic Acids Research* 32, e178. |
| `vermeulen1`, `vermeulen2` | Vermeulen, J. et al. (2009). External oligonucleotide standards enable cross laboratory comparison and exchange of real-time quantitative PCR data. *Nucleic Acids Research* 37, e138. |

For the remaining datasets (`competimer`, `dil4reps94`, `dyemelt`, `htPCR`, `reps`, `reps2`,
`reps3`, `reps384`, `testdat`), see the qpcR documentation.

If you use these datasets in published work, please cite qpcR and the original publication.
