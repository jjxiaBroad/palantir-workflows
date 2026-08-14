# Docker Setup

This directory contains the Docker image definition for `downsample.nf`'s processes (`SUMMARIZE_GEX_DOWNSAMPLE_TARGETS`, `SUMMARIZE_GRNA_DOWNSAMPLE_TARGETS`, `DOWNSAMPLE_MOLECULE_INFO`, `DOWNSAMPLE_CRISPR_ANNDATA`, `COMBINE_DOWNSAMPLED_BATCH`).

## Structure

```
docker/
├── README.md              # General Docker documentation
├── SETUP.md               # This file - build/push guide
├── build_and_push.sh      # Script to build and push the qc image to ECR
└── qc/
    ├── Dockerfile         # Single unified image for all downsample.nf processes
    └── requirements.txt   # Python package requirements (optional)
```

## Building and Pushing to ECR

### Prerequisites

1. **AWS CLI configured** with appropriate credentials:
   ```bash
   aws configure
   ```

2. **Docker installed** and running

3. **ECR Registry URL** - Find your ECR registry URL in the AWS Console or run:
   ```bash
   aws ecr describe-repositories --region us-east-1
   ```
   Format: `<account-id>.dkr.ecr.<region>.amazonaws.com`

### Build and Push

Set your ECR registry and build:

```bash
# Set your ECR registry URL
export ECR_REGISTRY=123456789012.dkr.ecr.us-east-1.amazonaws.com
export AWS_REGION=us-east-1

# Build and push with default 'latest' tag
./build_and_push.sh

# Or build with a specific tag (e.g., version or commit hash)
./build_and_push.sh v1.0.0
./build_and_push.sh $(git rev-parse --short HEAD)
```

The script will:
1. Authenticate Docker to your ECR registry
2. Create the ECR repository if it doesn't exist
3. Build the Docker image from `qc/Dockerfile`
4. Push the image to ECR

### Use the pushed image

`qc_container` is a required pipeline parameter (there is no config-file default) — pass it explicitly on the command line, or set `params.qc_container` in a config file included via `-c`:

```bash
nextflow run downsample.nf \
  --qc_container 123456789012.dkr.ecr.us-east-1.amazonaws.com/pipseq-downsample-qc:latest \
  --samplesheet samplesheet.csv \
  --batch_id "Batch_A" \
  --batch_basename "batch_a" \
  --outdir results \
  ...
```

## QC Image Contents

The `qc` image includes:
- **Python 3.10** with conda
- **Python packages**: pandas, numpy, scipy, h5py, anndata, scanpy (all via conda)
- **System tools**: build-essential

This single image is used by all processes in `downsample.nf`.

## Troubleshooting

**Authentication errors**: Ensure your AWS credentials are configured and have ECR push permissions

**Repository not found**: The script creates the repository automatically, but ensure your IAM role has `ecr:CreateRepository` permission

**Build failures**: Check that the Dockerfile path is correct and all dependencies are available

**Container fails to start / wrong container used**: Double check the `--qc_container` value passed on the command line (or set in your config file) — the pipeline has no default container image.
