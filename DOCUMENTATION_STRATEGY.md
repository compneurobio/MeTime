# MeTime Documentation Strategy (GitHub + pkgdown)

This document is a practical blueprint for making MeTime easier to understand, faster to onboard, and more compelling to adopt.

## 1) Documentation goals

MeTime docs should help three audiences:

1. **First-time users**: Understand what problem MeTime solves and run a first successful analysis in <10 minutes.
2. **Applied researchers**: Learn method-specific workflows (imputation, networks, regressions, etc.) with confidence.
3. **Contributors/developers**: Extend analysis modules without breaking class conventions.

A strong documentation system should:

- Reduce time-to-first-result.
- Make each analysis path discoverable and reproducible.
- Explain *when to use what* (method selection guidance).
- Increase trust with transparent assumptions, limitations, and validation examples.

---

## 2) Information architecture to adopt

Use this high-level structure in the repository:

```text
MeTime/
├── README.md
├── NEWS.md
├── CONTRIBUTING.md
├── CODE_OF_CONDUCT.md
├── docs/
│   ├── documentation-strategy.md
│   ├── getting-started/
│   │   ├── installation.md
│   │   ├── quickstart.md
│   │   └── sample-data-walkthrough.md
│   ├── user-guides/
│   │   ├── data-preparation.md
│   │   ├── building-pipelines.md
│   │   ├── choosing-methods.md
│   │   └── troubleshooting.md
│   ├── reference/
│   │   ├── function-families.md
│   │   ├── metime-analyser-class.md
│   │   └── plotting-system.md
│   ├── case-studies/
│   │   ├── disease-mechanism-example.md
│   │   └── multi-cohort-meta-analysis.md
│   └── developer/
│       ├── extending-modules.md
│       ├── testing-and-validation.md
│       └── release-checklist.md
└── vignettes/
    └── ...existing Rmd vignettes
```

### Why this helps

- `README.md` stays short and conversion-oriented.
- `docs/getting-started/*` gives a linear onboarding path.
- `docs/user-guides/*` answers practical “how do I …?” tasks.
- `docs/reference/*` explains architecture and naming conventions.
- `docs/case-studies/*` provides proof and inspiration.
- `docs/developer/*` lowers contribution friction.

---

## 3) README redesign (high-impact)

Your current README is rich but dense. Reframe it to be scan-friendly and action-first:

### Recommended README sections

1. **One-line value proposition**
   - “MeTime is an R package for modular longitudinal metabolomics analysis pipelines.”
2. **Who MeTime is for**
   - 2–3 bullets (biostatisticians, metabolomics researchers, translational teams).
3. **Why MeTime vs ad-hoc scripts**
   - Reproducibility, modularity, multi-dataset support, unified class.
4. **Install**
   - CRAN/GitHub install snippets.
5. **Quickstart (copy-paste runnable)**
   - Minimal pipeline using included sample data.
6. **Method map**
   - Table: method, question answered, key function(s), output type.
7. **Learning path**
   - “Start here” links to 3 foundational guides.
8. **Citation + support**
   - How to cite and where to ask questions.

### Conversion-oriented README elements

- Add badges: build status, coverage, CRAN/GitHub version, docs site.
- Include a single compact workflow image with alt text.
- Add “Common use cases” cards linking to vignettes.
- End with “What’s next” CTA: quickstart and first advanced guide.

---

## 4) Getting-started flow (first 10 minutes)

Create a dedicated onboarding path with strict success criteria:

### `docs/getting-started/quickstart.md`

Should include:

- Prerequisites and required R version.
- Install command and dependency notes.
- Load sample data.
- Build a minimal `metime_analyser` object.
- Run one `mod_*` + one `calc_*` + one plotting call.
- Show expected output snippets and one plot.
- “If this failed” section with top 5 errors.

### Success metric

A new user should reach first interpretable output in **<10 minutes** without reading deeper docs.

---

## 5) Method-oriented user guides

For each method family (distribution, feature selection, regression, networks, etc.), standardize vignette structure:

1. **When to use this method**
2. **Biological/statistical question answered**
3. **Input requirements**
4. **Pipeline steps**
5. **How to interpret results**
6. **Common pitfalls**
7. **Validation checks**
8. **Related methods / next step**

This avoids the common problem of users seeing code but not understanding decision criteria.

---

## 6) Function discoverability improvements

Users often struggle with many functions. Add a function-family index at `docs/reference/function-families.md`:

- `add_*`: enrich object with additional information.
- `mod_*`: transform/filter data.
- `calc_*`: run analysis modules.
- `meta_*`: compare/combine results.
- `get_*`: inspect/extract information.
- `write_*`: export outputs/reports.

For each family, provide:

- Intent
- Typical order of use
- 3–5 common functions
- Anti-patterns (what not to do)

---

## 7) Encourage usage with proof-oriented content

To motivate adoption, add social proof and scientific utility:

- **Case studies** with narrative outcomes (not just code).
- **Comparison table**: MeTime pipeline vs custom script workflow.
- **Reproducibility checklist** users can report in manuscripts.
- **Performance notes** (runtime/memory expectations).
- **“What questions can MeTime answer?”** matrix tied to study goals.

A user is more likely to adopt when they see direct mapping from package features to their research questions.

---

## 8) pkgdown site structure

Use `_pkgdown.yml` to mirror the architecture:

- Get started
- User guides
- Method tutorials
- Reference (grouped by function families)
- Developer docs
- Changelog/news

Also:

- Add clear navbar labels (“Get started”, “Methods”, “Reference”, “Developer”).
- Surface 3 top tutorials directly on homepage.
- Add a concise package lifecycle/status notice.

---

## 9) Quality standards for every doc page

Define a reusable template and enforce these checks:

- Plain-language summary in first 3 lines.
- Minimal runnable code chunk.
- Expected output shown.
- Input schema/assumption list.
- Troubleshooting section.
- Links to preceding and next pages.
- Last-updated date and package version compatibility.

---

## 10) Suggested 6-week execution roadmap

### Week 1: Foundation

- Rewrite README structure.
- Add `NEWS.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`.
- Create docs folders and page stubs.

### Week 2: Onboarding

- Build installation + quickstart + sample walkthrough pages.
- Add first-time-user troubleshooting.

### Week 3–4: Method guides

- Standardize 3 highest-usage vignettes first.
- Add method-selection guide (“choosing-methods”).

### Week 5: Reference + developer docs

- Build function-family index.
- Add module extension and testing docs.

### Week 6: polish + promotion

- Tighten internal links.
- Improve pkgdown navigation.
- Add case studies and usage encouragement section.

---

## 11) Documentation KPIs to track

Track impact like a product:

- README-to-quickstart click-through rate.
- Quickstart completion rate (if telemetry possible via examples/issues feedback).
- Number of “how do I start?” issues over time.
- Time-to-first successful pipeline (via user feedback form).
- Number of users referencing case studies in issues/discussions.

---

## 12) Immediate next steps for this repository

1. Keep the current detailed conceptual content, but split it across topical docs pages.
2. Convert README into an onboarding landing page (shorter + runnable).
3. Add a method map table and links to existing vignettes.
4. Create troubleshooting and function-family index pages before adding new features.
5. Ensure each vignette states assumptions and interpretation guidance clearly.

