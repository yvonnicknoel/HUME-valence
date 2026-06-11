# HUME-valence

Data, model code, and figure-reproduction script accompanying the paper:

> Noel, Y. (in press). *Emotional valence as a control variable in a generic valence-arousal-dominance dynamical system.* **Affective Science**, special issue *"The nature of valence: Pluses, minuses and progress"*.

The manuscript argues that emotional **valence** acts as a dynamic stress signal — a control
variable — rather than a static hedonic dimension, with arousal and dominance behaving as
nonlinear (single-peaked) functions of valence. Rather than standard factor analysis, the
psychometric approach relies on an **unfolding model** to recover a unidimensional valence
representation from emotion-induction data.

## Contents

| File | Description |
|------|-------------|
| `emotions-2024-2026.csv` | Experimental dataset (see below). |
| `BUM.R` | Implementation of the **Beta Unfolding Model** (BUM) as an `R6` class, with MAP/EM estimation. |
| `HUME_analyses_and_plots.R` | Analysis pipeline that reproduces Figures 5–9 of the paper. |

## The dataset

`emotions-2024-2026.csv` contains **326 participants** who rated **40 emotion items** on a
0–100 scale, **before** and **after** an emotion induction. Each item appears twice: a `1`
suffix for the pre-induction rating and a `2` suffix for the post-induction rating
(e.g. `Content1` / `Content2`).

Columns (83 total):

- `Gender`, `Age` — participant characteristics.
- `Induction` — induction condition, one of `Angry`, `Fear`, `Happy`, `Neutral`, `Sad`.
- `Content1` … `Aggressive1` — the 40 **pre-induction** emotion ratings.
- `Content2` … `Aggressive2` — the same 40 **post-induction** emotion ratings.

Induction-group sizes: Angry 62, Fear 72, Happy 67, Neutral 61, Sad 64.

## The model

`BUM.R` implements the Beta Unfolding Model introduced in:

> Noel, Y. (2014). A beta unfolding model for continuous bounded responses. *Psychometrika*, 79(4), 647–674.

This version adds **Gaussian priors** on the item parameters (`delta`, `lambda`, `tau`) for
MAP estimation, regularizing extreme item locations far from the person distribution.

Each item is characterized by a peak **location** (`delta`) along the valence dimension, an
**acceptability/log-slope** parameter (`lambda`), and a **dispersion** parameter (`tau`).

## Reproducing the figures

```r
# From within this directory, in an R session:
source("HUME_analyses_and_plots.R")
```

The script:

1. Loads `emotions-2024-2026.csv` and splits it into pre/post item blocks.
2. Sources `BUM.R` and fits the unfolding model to the pre-induction data, then to each
   post-induction condition (separately, then with item locations fixed to the pre-induction
   estimates so that person scores are comparable across time and condition).
3. Writes the figures as PDFs into an `Images/` subdirectory:

   | Figure | Output file | Content |
   |--------|-------------|---------|
   | Fig. 5 | `Images/Fig.5-BUM-params.pdf` | Illustration of the BUM item parameters. |
   | Fig. 6 | `Images/Fig.6-emotion15-BUM-pre.pdf` | 15 illustrative item characteristic curves (pre-induction). |
   | Fig. 7 | `Images/Fig.7-item-locations.pdf` | Emotion peak locations along the valence dimension, by domain. |
   | Fig. 8 | `Images/Fig.8-induction-effect.pdf` | Effect of each induction on participants' valence scores. |
   | Fig. 9 | `Images/Fig.9-ICC-impact.pdf` | Reactivity of response curves to the induction. |

   It also writes `subjects.csv` (participant pre/post valence scores).

> **Before running**, create the output directory the script writes to:
>
> ```r
> dir.create("Images", showWarnings = FALSE)
> ```
>
> The script also opens interactive `x11()` graphics windows; on systems without X11 you
> may wish to replace those calls with a suitable device.

### Dependencies

R (≥ 4.0) and the following CRAN packages:

```r
install.packages(c("R6", "tensor", "wordcloud", "Hmisc", "mgcv", "dglm"))
```

## Citation

If you use this dataset or code, please cite the *Affective Science* paper above and, for the
model itself, Noel (2014, *Psychometrika*).

## License

Released under the [MIT License](LICENSE).
