#!/usr/bin/env python
# -*- coding: utf-8 -*-
# import packages
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysam import VariantFile

# from snakemakefile the input, output, wildcards are attributes of snakemake
# object
quals = [record.qual for record in VariantFile(snakemake.input[0])]
plt.hist(quals)

plt.savefig(snakemake.output[0])