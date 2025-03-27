#!/bin/sh

# Unlock Snakemake in case of incomplete or failed runs
snakemake -s main.snk --unlock

# Run Snakemake with 5 jobs in parallel, each using 32 cores
snakemake -s main.snk all --cores 32 --jobs 1 --rerun-incomplete
