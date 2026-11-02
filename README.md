# Every Prime Is Born a Twin — publication package

Author: Anthony Phillips · 2 November 2026

| file | purpose | where it goes |
| --- | --- | --- |
| paper.pdf, paper.tex | research note: "The variance of twin primes across residue classes: a conditional asymptotic and numerical evidence" | Zenodo (upload PDF, get DOI), arXiv math.NT (upload .tex + fig*.pdf) |
| note.md, note.html | the same note as Markdown / standalone HTML (MathJax from CDN) | your site, GitHub |
| article.md, article.html | the expository article "Every Prime Is Born a Twin" with figures | your blog |
| fig1_brun.*, fig2_gpy.*, fig3_rigidity.* | figures (PDF for LaTeX, PNG for web) | both |
| figures.py | regenerates the figures from the numbers | — |
| scripts/ | all computations; hl2.c is the fast predictor (`gcc -O2 hl2.c -lm && ./hl2 1e9`) | GitHub |

Publishing order that gives you a defensible date: Zenodo first (DOI, minutes), then your site citing the DOI, then arXiv when you have an endorser, then the emails.

## Citation

Anthony Phillips, *The variance of twin primes across residue classes: a conditional asymptotic and numerical evidence*, preprint, 2 November 2026.
Zenodo DOI: 10.5281/zenodo.23891109  —  https://doi.org/10.5281/zenodo.23891109
Repository: https://github.com/justamancr/twin-primes

## Licenses

Code: MIT (`LICENSE`). Paper, article and figures: CC BY 4.0 (`LICENSE-CONTENT.md`).
