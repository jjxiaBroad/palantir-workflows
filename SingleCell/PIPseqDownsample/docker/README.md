# Docker and Container Management for Nextflow Pipeline

## Overview

Similar to WDL's `runtime.docker` attribute, Nextflow allows you to specify Docker containers at the **process level** via the `container` directive.

This pipeline uses one container image:
- **`qc_container`** (built from `docker/qc/Dockerfile` in this directory): used by every process in `downsample.nf` (`SUMMARIZE_GEX_DOWNSAMPLE_TARGETS`, `SUMMARIZE_GRNA_DOWNSAMPLE_TARGETS`, `DOWNSAMPLE_MOLECULE_INFO`, `DOWNSAMPLE_CRISPR_ANNDATA`, `COMBINE_DOWNSAMPLED_BATCH`). Bundles pandas/numpy/scipy/h5py/anndata/scanpy — no CRISPAT, Node, or other QC-report-only dependencies, since this pipeline doesn't run guide assignment or generate the HTML/Sankey QC reports that `../PIPseqPipeline` does.

It is a **required** pipeline parameter with no default — you must pass `--qc_container` explicitly.

## Directory Structure

```
docker/
├── README.md              # This file
├── SETUP.md               # Build/push guide for the qc image
├── build_and_push.sh      # Script to build and push the qc image to ECR
└── qc/
    ├── Dockerfile         # qc_container image definition
    └── requirements.txt   # Additional pip-only Python dependencies
```

## Quick Comparison: WDL vs Nextflow

**WDL:**
```wdl
task my_task {
    runtime {
        docker: "my-image:latest"
    }
}
```

**Nextflow:**
```groovy
process MY_PROCESS {
    container "my-image:latest"
}
```

See `SETUP.md` for how to build and push this pipeline's `qc_container` image to ECR.
