# Track A roadmap: tissue dynamics from a spatial snapshot

One of six connected tracks; see [scope-map.md](scope-map.md) for the others and how they link.

**Thesis.** How much tissue dynamics can one spatial snapshot actually identify? This
project reproduces OSDR, a 2025 *Nature* method that infers cell-population dynamics from a
single biopsy using Ki67 and neighbourhood composition. It then tests the method against
synthetic tissues with known dynamics and replaces cell-level statistics with patient-level
inference.

The goal is not a better OSDR. It is a clear, reproducible answer to when OSDR's conclusions
can be trusted.

Built from: Somer, Mannor & Alon, *Nature* 2025 (doi:10.1038/s41586-025-09876-1), and the
papers in `Sc-spatial-temporal-multi-omics/`. Each design choice below cites the paper it
comes from.

## What OSDR does, and where it is weak

OSDR models each cell type's population as

    dX_i/dt = X_i * ( p_i_plus(N) - p_i_minus )

`p_i_plus(N)` is the probability that a cell of type *i* is dividing. It is fitted as a
logistic regression of thresholded Ki67 on the counts of each cell type within 80 µm.
`p_i_minus` is the removal (death) rate. It is **not measured**: it is set to the mean
division rate, which is exact only at steady state. Fixed points and phase portraits follow
from the mean-field ODE.

These are the weaknesses this project tests. Each one maps to a CLAUDE.md rule.

| Weakness in the published analysis | Rule | What this project does |
|---|---|---|
| Death rate assumed constant, equal to mean division rate | Fail loudly; identifiability first | States that neighbourhood-dependent death is structurally unidentifiable from Ki67 alone. Quantifies the fixed-point bias with synthetic tissues where death does depend on neighbourhood |
| Statistics use cells as n (P < 10⁻¹³) | Replicates; the unit of replication is the patient | Patient-level bootstrap and mixed effects (patient and region intercepts). States the smallest detectable difference at ~50 patients per group |
| Responder and non-responder models fitted on outcome-defined groups | Baselines; no tuning against the answer | Leave-one-patient-out classification, reported with MCC and a confusion matrix |
| Radius 80 µm and Ki67 threshold fixed by hand | Never tune a threshold against the answer | Radius derived from Ki67-residual correlograms (Voyager); sensitivity scans over radius, threshold and neighbourhood type |
| The package plots a phase portrait even when assumptions are violated | Fail loudly | Refuses to draw a fixed point when state-space coverage or events-per-parameter is too low |
| All cohorts are breast-cancer IMC from overlapping groups | Generalize across ≥3 systems | Two open cohorts plus synthetic tissues, analysed leave-one-cohort-out |

## Systems

| System | Source | Access | Role |
|---|---|---|---|
| Synthetic tissues | This repository | n/a | Positive control and ceiling: known division and death rules, so recovery can be measured |
| Danenberg 2022, METABRIC IMC | Zenodo 7324285 | **Open**, CC BY 4.0, ~0.9 GB processed | Primary cohort; 715 patients; survival data |
| Fischer 2023 IMC | Zenodo 7494509 | **Open**, CC BY 4.0, ~39 GB | Second cohort, held out for replication |
| Wang 2023 NeoTRIP TNBC | Zenodo 7990870 | **Restricted** (application, academic non-commercial, 9.1 TB) | Optional. Never stored in or published from this repository |
| TCGA-BRCA bulk RNA-seq | NCI GDC | Open (processed) | Translation: deconvolve inferred tissue states, test against outcome |

Large data lives outside OneDrive (e.g. `C:\data\spatiotemporal`), set by an environment
variable.

## Milestones

Each milestone is a GitHub release with one figure and one table, regenerated from raw
inputs by one command.

### R0. Skeleton
`environment.yml` (conda-forge, every pin explained), `pyproject.toml`, ruff and pytest
via pre-commit, GitHub Actions running on a tiny synthetic tissue, `LICENSE`,
`CITATION.cff`.

### R1. Synthetic tissue simulator (built before any real data)
A 2D agent-based tissue whose dynamics are known exactly:
- Division and death rules as functions of neighbourhood, including death rules that
  depend on neighbourhood (which OSDR cannot see)
- Ki67 as a noisy observation with a finite "visible" window, not the division event itself
- Segmentation and phenotyping errors, worst for small immune cells (Deng, *Nat Rev
  Bioeng* 2023; Abdelaal, *Genome Biol* 2019: label errors flow straight into neighbour
  counts)
- About three regions per patient, with patient-level random effects

Release figure: a known phase portrait next to the one recovered from the simulator's snapshot.

### R2. Data layer
- Cached, checksummed fetch of the Danenberg cohort
- One validated cell table: patient, image, region, x and y in µm, phenotype, raw Ki67
  counts. Units and pixel size are read from the source, never inferred
- Show which images belong to which patient. Multiple cores share one patient ID and must not
  be counted as independent

### R3. Neighbourhoods
- Radius, k-nearest-neighbour and Delaunay neighbourhoods, with edge correction
- Radius derived from the spatial autocorrelation of Ki67 residuals rather than fixed
- Scan over neighbourhood definition (Gulati, *Nat Rev MCB* 2024)

### R4. Division model, with identifiability first
- Before fitting: events per parameter, a map of how well the data covers state space, and
  a profile over the unknown death scale
- Ki67 below the noise floor treated as censored rather than hard-thresholded
- Baselines it must beat on held-out patients: a constant rate per cell type, and Ki67
  labels permuted within each image and cell type (CPA, *Mol Syst Biol* 2023, found a
  linear baseline matched a deep model)
- Mixed-effects fit with patient and region intercepts

### R5. Dynamics
- Mean-field ODE, fixed points and Jacobian stability
- Phase portraits with patient-bootstrap confidence regions
- An explicit invalid result when coverage or identifiability checks fail

### R6. Reproduce one headline result
- The fibroblast–macrophage "hot" and "cold" fibrosis fixed points on Danenberg, compared
  against the official OSDR package as a reference implementation
- Replicate on Fischer, held out. Results are graded "robust" (holds in both) or
  "indication" (one cohort only), following the leave-one-dataset-out design of Muus,
  *Nat Med* 2021

### R7. Stress test (the new result)
On synthetic tissues, measure how far the recovered fixed points move when:
- death depends on neighbourhood,
- the tissue is not at steady state,
- segmentation error rises,
- tissue area per patient shrinks.

Report the ceiling as recovery using true division labels instead of Ki67. This answers
"when can OSDR be trusted?" with numbers.

### R8. Translation to cohorts
- Turn inferred tissue states into an expression signature, deconvolve TCGA-BRCA bulk
  data and test against survival with patient-level statistics. This follows Chaffin,
  *Nature* 2022, which validated a single-nucleus cell state by deconvolving independent
  bulk cohorts
- Re-test OSDR's survival association (log-rank P = 0.0046, n = 607) against simple
  baselines: tumour fraction and mean Ki67

### R9. Wet-lab layer and experiment planner
- Per-cohort notebooks following `docs/notebook-design.md`: question card, status lines, a
  bench glossary (division rate per day, the Ki67 window), and baseline next to every fit
- **Design your experiment.** Given a planned IMC study, how much tissue area and how many
  patients are needed to identify a fixed point, and which added marker narrows the
  uncertainty most? The expected answer *was* a death or apoptosis marker (e.g. cleaved
  caspase-3), which would make the removal rate measurable instead of assumed. **Tested
  2026-10-28 and it did not hold:** in the Rendeiro COVID-19 lung data, cleaved caspase-3
  is higher in macrophages (0.62) and neutrophils (0.57) than in the authors' own "Dying
  cells" cluster (0.40), and gives an interior fixed point in 0 of 23 patients. A panel
  slot is not a measurement. The planner must therefore rank a *candidate* marker by
  whether it separates the population it names in that tissue, and report that separation
  next to the expected information gain, rather than assuming it. See the README's
  negative results. The planner ranks candidate experiments by expected information gain.
  An optional language-model step drafts the protocol rationale

### R10. Packaging
Snakemake workflow, documentation site, Zenodo DOI.

## How public databases fit

| Database | Use | Counts as validation? |
|---|---|---|
| TCGA (GDC), METABRIC | Bulk deconvolution and survival (R8) | Yes, independent cohorts |
| CELLxGENE, HCA | Reference cell-type signatures for deconvolution | No, input |
| Human Protein Atlas, GTEx | Sanity check that IMC panel targets are expressed where expected | No, annotation |
| GWAS Catalog, gnomAD | Breast-cancer risk loci or variant annotation, for prioritization only | No. None of the papers in this collection used them to validate dynamics |

## Role coverage

| Role | Covered by |
|---|---|
| MD Anderson, Immunotherapy Platform data scientist | Spatial single-cell analysis of the tumour microenvironment, cell–cell interaction dynamics, TCGA, reusable GitHub pipelines (R2–R8) |
| BMS, Translational Computational Biology | Spatial proteomics, a response biomarker evaluated with patient-level validation, deconvolution, reproducible workflows (R4–R8, R10). The posting's PhD + 4 years requirement is not something a project can substitute |
| LBNL, Quantitative Modeling | A model that proposes the next experiment (R9) |

## Links to other tracks

- **Track B (infection inside one host).** Uses the same ODE, identifiability and fixed-point
  code. Infected-cell dynamics and tissue cell dynamics are the same mathematics.
- **Track E (variant to cell type).** Asks whether breast-cancer GWAS risk-gene programs
  concentrate in cell types whose division depends on neighbourhood. This is annotation, and is
  labelled as such.
- **Track F (experiment planner).** R9's planner is the Track A instance of the shared planner.
- **Methods borrowed from the cardiac and SARS-CoV-2 papers.** Donor-level pseudo-bulk and
  leave-one-dataset-out testing (Muus 2021); bulk-deconvolution validation (Chaffin 2022).

## Not in this track

- Training on or redistributing anything derived from the restricted NeoTRIP data
