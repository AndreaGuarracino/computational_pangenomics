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

The [human leukocyte antigen (HLA)](https://en.wikipedia.org/wiki/Human_leukocyte_antigen) system is a complex of genes on chromosome 6 in humans which encode cell-surface proteins responsible for the regulation of the immune system.

Let's build a pangenome graph from a collection of sequences of the DRB1-3123 gene downloaded from [HLA-zoo](https://github.com/ekg/HLA-zoo):

```shell
pggb -i HLA-zoo/seqs/DRB1-3123.fa -n 12 -t 8 -o DRB1_3123.1
```

Run `pggb` without parameters to get information on the meaning of each parameter:

```shell
pggb

usage: pggb -i <input-fasta> -o <output-dir> [options]
options:
   [wfmash]
   -i, --input-fasta FILE      input FASTA/FASTQ file
   -s, --segment-length N      segment length for mapping [default: 5000]
   -l, --block-length N        minimum block length filter for mapping [default: 5*segment-length]
   -p, --map-pct-id PCT        percent identity for mapping/alignment [default: 90]
   -c, --n-mappings N          number of mappings for each segment [default: 1]
   -g, --hg-filter-ani-diff N  filter out mappings unlikely to be N% less than the best mapping [default: %]
   -N, --no-splits             disable splitting of input sequences during mapping [default: enabled]
   -x, --sparse-map N          keep this fraction of mappings ('auto' for giant component heuristic) [default: 1.0]
   -K, --mash-kmer N           kmer size for mapping [default: 19]
   -F, --mash-kmer-thres N     ignore the top % most-frequent kmers [default: 0.001]
   -Y, --exclude-delim C       skip mappings when the query and target have the same
                              prefix before the last occurrence of the given character C [default: #, assuming PanSN-spec]
   [seqwish]
   -k, --min-match-len N       filter exact matches below this length [default: 23]
   -f, --sparse-factor N       keep this randomly selected fraction of input matches [default: no sparsification]
   -B, --transclose-batch      number of bp to use for transitive closure batch [default: 10M]
   [smoothxg]
   -X, --skip-normalization    do not normalize the final graph [default: normalize the graph]
   -n, --n-haplotypes N        number of haplotypes
   -j, --path-jump-max         maximum path jump to include in block [default: 0]
   -e, --edge-jump-max N       maximum edge jump before breaking [default: 0]
   -G, --poa-length-target N,M target sequence length for POA, one per pass [default: 700,900,1100]
   -P, --poa-params PARAMS     score parameters for POA in the form of match,mismatch,gap1,ext1,gap2,ext2
                              may also be given as presets: asm5, asm10, asm15, asm20
                              [default: 1,19,39,3,81,1 = asm5]
   -O, --poa-padding N         pad each end of each sequence in POA with N*(mean_seq_len) bp [default: 0.001]
   -d, --pad-max-depth N       depth/haplotype at which we don't pad the POA problem [default: 100]
   -b, --run-abpoa             run abPOA [default: SPOA]
   -z, --global-poa            run the POA in global mode [default: local mode]
   -M, --write-maf             write MAF output representing merged POA blocks [default: off]
   -Q, --consensus-prefix P    use this prefix for consensus path names [default: Consensus_]
   [odgi]
   -v, --skip-viz              don't render visualizations of the graph in 1D and 2D [default: make them]
   -S, --stats                 generate statistics of the seqwish and smoothxg graph [default: off]
   [vg]
   -V, --vcf-spec SPEC         specify a set of VCFs to produce with SPEC = REF[:LEN][,REF[:LEN]]*
                              the paths matching ^REF are used as a reference, while the sample haplotypes
                              are derived from path names, assuming they match the PanSN; e.g. '-V chm13',
                              a path named HG002#1#ctg would be assigned to sample HG002 phase 1.
                              If LEN is specified and greater than 0, the VCFs are decomposed, filtering 
                              sites whose max allele length is greater than LEN. [default: off]
   [multiqc]
   -m, --multiqc               generate MultiQC report of graphs' statistics and visualizations,
                              automatically runs odgi stats [default: off]
   [general]
   -o, --output-dir PATH       output directory
   -D, --temp-dir PATH         directory for temporary files
   -a, --input-paf FILE        input PAF file; the wfmash alignment step is skipped
   -r, --resume                do not overwrite existing outputs in the given directory
                              [default: start pipeline from scratch]
   -t, --threads N             number of compute threads to use in parallel steps [default: 16]
   -T, --poa-threads N         number of compute threads to use during POA (set lower if you OOM during smoothing)
   -A, --keep-temp-files       keep intermediate graphs
   -Z, --compress              compress alignment (.paf), graph (.gfa, .og), and MSA (.maf) outputs with pigz,
                              and variant (.vcf) outputs with bgzip
   --names-with-params         put parameter values in filenames, instead of hashes
   --version                   display the version of pggb
   -h, --help                  this text

Use wfmash, seqwish, smoothxg, odgi, gfaffix, and vg to build, project and display a pangenome graph.
```

Take a look at the files in the `DRB1_3123.1` folder.

We get a graph in GFA (`*.gfa`) and ODGI (`*.og`) formats. These can be used downstream in many tools, including those in `vg`. You can visualize the GFA format graph with [`BandageNG`](https://github.com/asl/BandageNG), and use `odgi` directly on the `*.gfa` or `*.og` output.

### Understanding `odgi` visualizations

We obtain a series of diagnostic images that represent the pangenome alignment. These are created with `odgi viz` (1D matrix) and `odgi layout` with `odgi draw` (2D graph drawings).

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
