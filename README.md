## Getting started

```
snakemake -j 1 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=100000 --cluster "sbatch -t 10080 -J camip -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k -n
```
