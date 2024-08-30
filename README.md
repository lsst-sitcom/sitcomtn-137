[![Website](https://img.shields.io/badge/sitcomtn--137-lsst.io-brightgreen.svg)](https://sitcomtn-137.lsst.io)
[![CI](https://github.com/lsst-sitcom/sitcomtn-137/actions/workflows/ci.yaml/badge.svg)](https://github.com/lsst-sitcom/sitcomtn-137/actions/workflows/ci.yaml)

# Getting Started with Cell-Based Coadds

## SITCOMTN-137

As development for cell-based coadds continue, their testing will become pertinent during the commissioning process. This technote is meant to be an initial guide to generating and using cell-based coadds within the context of the LSST Science Pipelines and USDF, where several example use cases will be outlined in the form of brief code snippets and  initial analyses. Example code will also be maintained in the form of Jupyter notebooks, currently found [here](https://github.com/mirarenee/notebooks/tree/main/cell_coadds/technote).

**Links:**

- Publication URL: https://sitcomtn-137.lsst.io
- Alternative editions: https://sitcomtn-137.lsst.io/v
- GitHub repository: https://github.com/lsst-sitcom/sitcomtn-137
- Build system: https://github.com/lsst-sitcom/sitcomtn-137/actions/


## Build this technical note

You can clone this repository and build the technote locally if your system has Python 3.11 or later:

```sh
git clone https://github.com/lsst-sitcom/sitcomtn-137
cd sitcomtn-137
make init
make html
```

Repeat the `make html` command to rebuild the technote after making changes.
If you need to delete any intermediate files for a clean build, run `make clean`.

The built technote is located at `_build/html/index.html`.

## Publishing changes to the web

This technote is published to https://sitcomtn-137.lsst.io whenever you push changes to the `main` branch on GitHub.
When you push changes to a another branch, a preview of the technote is published to https://sitcomtn-137.lsst.io/v.

## Editing this technical note

The main content of this technote is in `index.md` (a Markdown file parsed as [CommonMark/MyST](https://myst-parser.readthedocs.io/en/latest/index.html)).
Metadata and configuration is in the `technote.toml` file.
For guidance on creating content and information about specifying metadata and configuration, see the Documenteer documentation: https://documenteer.lsst.io/technotes.
