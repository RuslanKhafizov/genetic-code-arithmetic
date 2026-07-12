# Project Rules

## Project scope

This repository supports the iterative development of the preprint
"Arithmetic structure of the standard genetic code" and its reproducible
data and code. The long-term goal is a journal-quality submission, currently
targeting *BioSystems* or a comparable journal.

Scientific correctness, logical dependency, transparent null models, and
reproducibility take priority over cosmetic editing. Do not overstate the
biological interpretation of an arithmetic result.

## Language and authoritative documents

- Communicate with the project owner in Russian unless asked otherwise.
- The Russian manuscript is the leading version: `russian/main_ru.tex`.
- The Russian README is `russian/README_ru.md`.
- The English manuscript `main.tex` and English `README.md` are synchronized
  translations intended for international publication.
- For a substantive manuscript or README change, edit the Russian version
  first, then update the English version in the same branch. Preserve
  scientific meaning rather than translating mechanically.
- Code, CSV files, metadata, website files, commit messages, pull requests,
  and release notes remain English-only.

## Scientific levels and claim discipline

- Always state which level is being analyzed:
  1. the 60-codon sense pool, excluding the three stops and initiator ATG;
  2. Key 0;
  3. Key 1.
- Do not treat a Key 0 or Key 1 result as assignment-free.
- Distinguish explicitly among:
  - primitive empirical facts;
  - definitions or imposed assignments;
  - algebraic consequences;
  - repeated copies propagated through nested codon groups.
- When reporting a family of identities, provide its dependency structure or
  a minimal independent generating set whenever possible.
- Do not count dependent equalities or divisibilities as independent evidence.
- Label pre-specified tests as confirmatory and post-hoc additions as
  exploratory or sensitivity analyses.
- Do not infer a biological mechanism, evolutionary cause, intentional
  encoding, or uniqueness over all conceivable genetic codes without a
  separate argument and evidence.

## Data and code

- Treat scripts as the authoritative definitions of generated results.
- Do not hand-edit a generated CSV to make it agree with prose. Correct the
  input or generating script, regenerate the output, and document the change.
- Preserve exact integer arithmetic and deterministic ordering where the
  existing scripts use them.
- Keep scripts standard-library-only unless a new dependency is scientifically
  necessary and explicitly approved.
- Keep input data, generated data, deterministic controls, Monte Carlo output,
  and independent verification conceptually separate.

## Working method

- Start each coherent task from an up-to-date `main` branch and use one
  descriptive branch for that task.
- Inspect the current branch, working-tree status, and diff before editing.
- Preserve unrelated user changes; never discard or overwrite them silently.
- Keep scientific changes, mechanical file moves, and release housekeeping in
  separate commits or pull requests.
- Do not reorganize repository paths together with changes to scientific logic.
- Do not commit, push, merge, tag, publish a GitHub Release, or modify Zenodo
  unless the project owner explicitly authorizes that step.
- The project owner uses GitHub Desktop and compiles both PDFs in Overleaf with
  pdfLaTeX. Do not replace tracked PDFs unless newly compiled files have been
  supplied or the owner explicitly requests local compilation.
- Only one AI assistant may edit the working tree at a time. A second assistant
  may review the diff read-only. Finish or commit a task before handing editing
  control to another assistant.

## Verification commands

Run commands from the repository root with an available Python 3.8+ interpreter.
The scripts require only the Python standard library.

- Core deterministic outputs: `python reproduce.py`
- Deterministic controls: `python controls.py`
- Monte Carlo result of record: `python mc_significance.py`
- Independent Monte Carlo cross-check: `python verify_mc.py`

`mc_significance.py` uses 1,000,000 trials and seed 0 by default. Run the full
Monte Carlo and independent verification when their code, inputs, statistics,
or reported results change. Do not rerun expensive stochastic checks for a
documentation-only change when the underlying files are unchanged.

After a generator is run, inspect the complete diff and confirm that only the
expected generated files changed. For a release, also verify that the English
and Russian PDFs open correctly and contain the synchronized results.

## Release model

- The GitHub/Zenodo Software line contains the complete repository, including
  code, CSV files, TeX sources, both PDFs, and the `russian/` directory.
- Publishing a GitHub Release automatically creates the corresponding Zenodo
  Software version. Treat GitHub's **Publish release** action as the final
  publication step for both services.
- The separate Zenodo preprint line contains only the English PDF and has an
  independent version history. Do not update it unless explicitly requested.
- Before a release, check the version and date in `CITATION.cff`, both README
  files, both TeX sources, both PDFs, release notes, and the exact target commit.

## Handoff requirements

At the end of a task, report:

- the active branch;
- files changed;
- checks run and their outcomes;
- any generated files or PDFs still requiring user action;
- deferred scientific questions and known limitations.
