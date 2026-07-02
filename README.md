## HICAR ➡️ SCHNAPS

Since July 2026, HICAR development has shifted to the [SCHNAPS](https://codeberg.org/SCHNAPS-Model/SCHNAPS) model. 
SCHNAPS is the successor to HICAR, running on GPUs, being easier to build and setup runs, and featuring improved snow-atmosphere coupling.
Also, the name is better.

You can find SCHNAPS now hosted on CodeBerg here: [https://codeberg.org/SCHNAPS-Model/SCHNAPS](https://codeberg.org/SCHNAPS-Model/SCHNAPS).
The GitHub repo is now archived and will not be updated.

### Model Users (No code changes)

To migrate a local clone of the model, simply:

```bash
git remote set-url origin https://codeberg.org/SCHNAPS-Model/SCHNAPS.git
git fetch origin
```

If you don't have working changes in your clone, you may just consider cloning everything fresh.

### Model Developers (Possible code changes)

#### Starting from a new fork

1. Create a [Codeberg account](https://codeberg.org/user/sign_up) and fork the new
   repository at https://codeberg.org/SCHNAPS-Model/SCHNAPS.

2. Point your local clone: `origin` at your new Codeberg fork, `upstream` at the project:

```bash
   git remote set-url origin   https://codeberg.org/YOUR-USERNAME/SCHNAPS.git
   git remote set-url upstream https://codeberg.org/SCHNAPS-Model/SCHNAPS.git
   git fetch --all
```

#### Migrating your local clone's changes

Assuming you are on some branch, `my_branch`:

```bash
v1_commit=d4305d4d12

# 0. **If you have changes in the working tree:** Preserve the dirty tree as a WIP commit + a backup branch
git add -A # or a less aggressive, thoughtful add...
git commit -m "WIP: local uncommitted changes"
WIP=$(git rev-parse HEAD); git branch wip_backup
git reset --hard HEAD~1

# 1. Point remotes at Codeberg and fetch upstream
git remote set-url origin   https://codeberg.org/YOUR-USERNAME/SCHNAPS.git
git remote remove upstream 2>/dev/null
git remote add upstream https://codeberg.org/SCHNAPS-Model/SCHNAPS.git
git fetch upstream

# 2. Merge SCHNAPS v1.0.0 onto your branch.
git merge "$v1_commit"
git commit            # conclude the merge once conflicts are resolved

#3. Collapse to a clean base
git reset --soft "$v1_commit"
git commit -m "my_branch: work ported onto SCHNAPS main"

# 4. Verify your changes did not add old hicar naming to anything, then commit (skip the commit if nothing changed).
git grep -i hicar || echo "clean: no 'hicar' left in tracked content"
git ls-files | grep -i hicar || echo "clean: no 'hicar' left in tracked paths"
git commit -am "Port my_branch additions to SCHNAPS naming" || echo "nothing to rename"

# 5. **If you had changes in the working tree:** Restore the originally-uncommitted changes as UNSTAGED, rename-aware
git cherry-pick -n "$WIP"     # uses WIP's parent as base → isolates only the dirty diff
#    If that dirty diff referenced HICAR names, re-run step 4a before continuing.
git reset                     # unstage, so it matches the original dirty shape

# 6. Publish to your fork when happy
git push -u origin my_branch
```

# The High-resolution Intermediate Complexity Atmospheric Research Model (HICAR)

<!-- The four merge-gate lanes that sign main's HEAD: hicar-full-test, valgrind, gpu-check, snow-parity -->
[![full-test](https://github.com/HICAR-Model/HICAR/actions/workflows/hicar-full-test.yml/badge.svg?branch=main)](https://github.com/HICAR-Model/HICAR/actions/workflows/hicar-full-test.yml)
[![valgrind](https://github.com/HICAR-Model/HICAR/actions/workflows/valgrind-memcheck.yml/badge.svg?branch=main)](https://github.com/HICAR-Model/HICAR/actions/workflows/valgrind-memcheck.yml)
[![GPU](https://github.com/HICAR-Model/HICAR/actions/workflows/gpu.yml/badge.svg?branch=main)](https://github.com/HICAR-Model/HICAR/actions/workflows/gpu.yml)
[![SNOWPACK parity](https://github.com/HICAR-Model/HICAR/actions/workflows/snowpack-compare.yml/badge.svg?branch=main)](https://github.com/HICAR-Model/HICAR/actions/workflows/snowpack-compare.yml)

[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![version](https://img.shields.io/github/v/tag/HICAR-Model/HICAR?filter=v*&label=version)](https://github.com/HICAR-Model/HICAR/tags)
[![DOI](https://zenodo.org/badge/638935780.svg)](https://zenodo.org/badge/latestdoi/638935780)
[![Paper: GMD](https://img.shields.io/badge/paper-GMD-1f7a8c.svg)](https://doi.org/10.5194/gmd-2023-16)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue.svg)](docs/index.md)

HICAR is a variant of the Intermediate Complexity Atmospheric Research (ICAR) model developed for sub-kilometer resolutions. The model is developed for downscaling of kilometer-scale NWP model output to resolutions used for land-surface simulations. HICAR features physics parameterizations shared by traditional weather models such as WRF, but with massively simplified dynamics which enable run times up to roughly two to three orders of magnitude faster than WRF (Reynolds et al., 2023).

#### Reference

Reynolds, D. S., Gutmann, E., Kruyt, B., Haugeneder, M., Jonas, T., Gerber, F., Lehning, M., and Mott, R.: The High-resolution Intermediate Complexity Atmospheric Research (HICAR v1.0) Model Enables Fast Dynamic Downscaling to the Hectometer Scale, Geosci. Model Dev. Discuss. [preprint], https://doi.org/10.5194/gmd-2023-16, in review, 2023. 

Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen (2016), *The Intermediate Complexity Atmospheric Research Model*, J. Hydrometeor, doi:[10.1175/JHM-D-15-0155.1](http://dx.doi.org/10.1175/JHM-D-15-0155.1).
