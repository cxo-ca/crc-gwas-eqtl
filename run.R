# Public-data pipeline. Run from the project root.
# Use Rscript run.R --offline to require the archived API responses.
source('scripts/01_fetch_public.R')
source('scripts/02_analyze_public.R')
