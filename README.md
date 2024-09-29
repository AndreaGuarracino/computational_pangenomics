# Computational Pangenomics

## Overview

Welcome to the Computational Pangenomics repository! This two-day workshop, held on October 7-8, 2024, is designed to introduce participants to the hot field of pangenomics and provide hands-on experience with pangenome analysis tools.

## Instructor

**Andrea Guarracino, PhD**  
Department of Genetics, Genomics and Informatics,
University of Tennessee Health Science Center
Memphis, TN, USA

Curriculum vitae at https://andreaguarracino.github.io/

## Workshop Details

- **Date**: October 7-8, 2024
- **Location**: Virtual

For more information, visit the [official workshop page](https://www.afinsubria.org/2024/07/09/computational-pangenomics/).

## Program

### Day 1: October 7, 2024

1. **Lecture**: Introduction to Pangenomics
   - Basic concepts and terminology
   - Overview of pangenome graph structures

2. **Practical**: Building Pangenome Graphs
   - Step-by-step guide to building simple pangenome graphs
   - Hands-on exercise with real data

### Day 2: October 8, 2024

3. **Lecture**: Understanding Pangenomes
   - Advanced pangenome concepts
   - Applications in population genetics and evolution

4. **Practical**: Analyzing Pangenomes
   - Techniques for pangenome analysis
   - Visualization tools for pangenome graphs

## Warming-up

### Build HLA pangenome graphs

The [human leukocyte antigen (HLA)](https://en.wikipedia.org/wiki/Human_leukocyte_antigen) system is a complex of genes on chromosome 6 in humans that encode cell-surface proteins responsible for regulating the immune system.

Let's build a pangenome graph from a collection of sequences of the DRB1-3123 gene downloaded from [HLA-zoo](https://github.com/ekg/HLA-zoo):

```shell
pggb -i HLA-zoo/seqs/DRB1-3123.fa -n 12 -t 8 -o DRB1_3123.1
```

To get information on the meaning of each parameter, run `pggb` without parameters:

```shell
pggb

usage: pggb -i <input-fasta> -o <output-dir> [options]
options:
   [wfmash]
   -i, --input-fasta FILE      input FASTA/FASTQ file
   -s, --segment-length N      segment length for mapping [default: 5000]
   ...
# The full output has been omitted for brevity.
```

Examine the files in the `DRB1_3123.1` folder. You'll find a graph in GFA (`*.gfa`) and ODGI (`*.og`) formats. These can be used downstream in many tools, including those in [`vg`](https://github.com/vgteam/vg). You can visualize the GFA format graph with [`BandageNG`](https://github.com/asl/BandageNG), and use [`odgi`](https://github.com/pangenome/odgi) directly on the `*.gfa` or `*.og` output.

### Understanding `odgi` visualizations

`odgi` generates a series of diagnostic images that represent the pangenome alignment. These are created with `odgi viz` (1D matrix) and `odgi layout` with `odgi draw` (2D graph drawings).

First, the 2D layout gives us a view of the total alignment. For small graphs, we can look at the version that shows where specific paths go (`*.draw_multiqc.png`):

![draw_multiqc.png](DRB1_3123.1/DRB1-3123.fa.bf3285f.11fba48.9c6ea4f.smooth.final.og.lay.draw_multiqc.png)

For larger ones, the `*.draw.png` result is usually more legible, but it lacks path information:

![draw.png](DRB1_3123.1/DRB1-3123.fa.bf3285f.11fba48.9c6ea4f.smooth.final.og.lay.draw.png)

We also get some 1D visualizations. Across the x-axis we have nodes of the graph (scaled by length) and across the y-axis we have paths (the sequences) that have been embedded in the graph.

This layout is capable of representing several kinds of information using color.

The default associates a color with each path. This is stable across different runs of `odgi viz`:

![viz_multiqc.png](DRB1_3123.1/DRB1-3123.fa.bf3285f.11fba48.9c6ea4f.smooth.final.og.viz_multiqc.png)

We also have a view that shows the "self depth" across the graph.
In this case there are no looping paths, so the color is always gray=1x.

![viz_depth_multiqc.png](DRB1_3123.1/DRB1-3123.fa.bf3285f.11fba48.9c6ea4f.smooth.final.og.viz_depth_multiqc.png)

We can look at orientation of paths using two views.

One shows the "position" of each path relative to the graph. It runs light to dark from 0 to path length.

![viz_pos_multiqc.png](DRB1_3123.1/DRB1-3123.fa.bf3285f.11fba48.9c6ea4f.smooth.final.og.viz_pos_multiqc.png)

A similar view shows inverted regions of paths relative to the graph in red, while the forward orientation in black.

![viz_inv_multiqc.png](DRB1_3123.1/DRB1-3123.fa.bf3285f.11fba48.9c6ea4f.smooth.final.og.viz_inv_multiqc.png)

And finally, a compressed view shows coverage across the pangenome coordinate space of all paths. It's a kind of heatmap. This helps when we have a lot of paths to consider:

![viz_O_multiqc.png](DRB1_3123.1/DRB1-3123.fa.bf3285f.11fba48.9c6ea4f.smooth.final.og.viz_O_multiqc.png)

### Graph statistics

Use `odgi stats` to obtain the graph length, and the number of nodes, edges, and paths:

```shell
odgi stats -i DRB1_3123.1/DRB1-3123.fa.bf3285f.11fba48.9c6ea4f.smooth.final.og

#length	nodes	edges	paths	steps
22237	4735	6481	12	34047
```

To determine if the resulting pangenome graph represents the input sequences well, compare the graph length and number of paths to the length and number of the input sequences.

### The effect of the minimum match filter `-k`

The `-k` parameter affects the behavior of `seqwish`. This filter removes exact matches from alignments that are shorter than `-k`. Short matches occur in regions of high diversity. We remove them to simplify the base graph structure.

Let's try setting a much higher `-k` than the default (`-k 23`):

```shell
pggb -i HLA-zoo/seqs/DRB1-3123.fa -n 12 -k 47 -t 8 -o DRB1_3123.2
```

The graph starts to become underaligned:

![draw_multiqc.png](DRB1_3123.2/DRB1-3123.fa.bf3285f.7608fc1.9c6ea4f.smooth.final.og.lay.draw_multiqc.png)

![viz_multiqc.png](DRB1_3123.2/DRB1-3123.fa.bf3285f.7608fc1.9c6ea4f.smooth.final.og.viz_multiqc.png)

### The effect of minimum pairwise identity `-p` of homology mapping

The `-p` setting affects the level of pairwise divergence accepted in the mapping step. Let's set this higher than the default (`-p 90`):

```shell
pggb -i HLA-zoo/seqs/DRB1-3123.fa -p 95 -n 12 -t 8 -o DRB1_3123.3
```

We lose mappings and this is visible in the diagnostic plots, which show that the graph has been broken into isolated components formed by sets of sequences that have >95% pairwise identity:

![draw_multiqc.png](DRB1_3123.3/DRB1-3123.fa.35d2267.11fba48.3a8f1bc.smooth.final.og.lay.draw_multiqc.png)

![viz_multiqc.png](DRB1_3123.3/DRB1-3123.fa.35d2267.11fba48.3a8f1bc.smooth.final.og.viz_multiqc.png)

### The effect of mapping segment length `-s`

Pangenome variation graphs built by `pggb` are based on homology mappings built using segments of a fixed size, rather than short k-mers (as in `minimap2`). This makes them suitable for quickly finding high-level patterns of homology.
You can think of `-s` as a seed length for the mappings.
It defaults to `5kb`, which testing has shown to provide a good tradeoff for computational efficiency, graph collinearity, and structural variant breakpoint detection.
Setting `-s` much higher can start to reduce sensitivity to small homologies.
However, when running with large eukaryotic genomes, we often set `-s` higher, sometimes up to `50k`. This can make the graph construction much more tractable.

### Call variants from the pangenome

We can write the variants in the input pangenome to a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format file. For this, pggb uses the vg deconstruct command to project the graph (in GFA file format) into VCF format. To invoke the projection, we need to specify which sequence to use as the reference for calling variants. Let's use the first sequence as an example:


```shell
pggb -i HLA-zoo/seqs/DRB1-3123.fa -p 95 -n 12 -t 8 -o DRB1_3123.4 -V 'gi|568815592'

[pggb] warning: there are sequence names (like 'gi|568815592:32578768-32589835') that do not match the Pangenome Sequence Naming (PanSN).
[pggb] ERROR: -V/--vcf-spec cannot be used if the Pangenome Sequence Naming (PanSN) is not respected.
```

We receive an error because calling variants requires modifying the sequence names to follow the [PanSN](https://github.com/pangenome/PanSN-spec) convention.
Let's hack the file:

```shell
sed 's/:/#1#/g' HLA-zoo/seqs/DRB1-3123.fa > HLA-zoo/seqs/DRB1-3123.pansn.fa
samtools faidx HLA-zoo/seqs/DRB1-3123.pansn.fa # index the new FASTA file
```

Now, let's try again:

```shell
pggb -i HLA-zoo/seqs/DRB1-3123.pansn.fa -p 95 -n 12 -t 8 -o DRB1_3123.4 -V 'gi|568815592'
```

To count how many variants are present in the VCF file:

```shell
grep '^#' -v DRB1_3123.4/DRB1-3123.pansn.fa.35d2267.11fba48.3a8f1bc.smooth.final.gi_568815592.vcf -c

2
```

We've identified only two variants! This low count is due to using `-p 95`, which is not suitable for this highly diverse region. Let's use the default `-p 90` instead:

```shell
pggb -i HLA-zoo/seqs/DRB1-3123.pansn.fa -n 12 -t 8 -o DRB1_3123.5 -V 'gi|568815592'

grep '^#' -v DRB1_3123.5/DRB1-3123.pansn.fa.bf3285f.11fba48.9c6ea4f.smooth.final.gi_568815592.vcf -c

1157
```
