# syri2jbrowse

Look at [SyRI](https://github.com/schneebergerlab/syri)'s structural
rearrangements between genomes in [JBrowse 2](https://jbrowse.org): each
syntenic, inverted, translocated and duplicated region a ribbon colored by its
type, the genomes stacked the way
[plotsr](https://github.com/schneebergerlab/plotsr) stacks them, in a view you
can zoom, search and add your own tracks to.

It takes the files plotsr takes, so a plotsr command becomes a JBrowse config by
changing its name.

## Install

```bash
pip install git+https://github.com/cmdcolin/syri2jbrowse
```

It needs nothing but Python 3.9 or later.

## Run

The command you ran plotsr with:

```bash
plotsr --sr A_Bsyri.out --sr B_Csyri.out --genomes genomes.txt -o plot.png
```

runs as:

```bash
syri2jbrowse --sr A_Bsyri.out --sr B_Csyri.out --genomes genomes.txt -o syri_jbrowse
```

and writes `syri_jbrowse/config.json` beside the files it reads. Open that config
in [JBrowse Desktop](https://jbrowse.org/jb2/docs/quickstart_desktop/)
(**Open config.json or .jbrowse file...**, or
`jbrowse-desktop syri_jbrowse/config.json`), or serve the folder beside
[JBrowse Web](https://jbrowse.org/jb2/docs/quickstart_web/) and open
`?config=syri_jbrowse/config.json`. It opens on every genome stacked, regions
under 10 kb hidden as plotsr hides them (`-s`, the same flag).

### Every genome in the first one's coordinates

plotsr compares each genome with the one above it. Give one more SyRI run per
genome against the first, and the config also draws every genome as a lane in
the first genome's coordinates, with a row of its SyRI regions on that genome:

```bash
syri2jbrowse --genomes genomes.txt \
  --sr A_Bsyri.out --sr B_Csyri.out --sr C_Dsyri.out \
  --pair A C A_Csyri.out --pair A D A_Dsyri.out
```

In a linear genome view of A, the track **SyRI regions** opens as lanes and
**SyRI regions on A, by genome** as a row per genome.

## Inputs

- `--genomes`: plotsr's genomes file, one `<file>\t<name>[\t<tags>]` line per
  genome in stack order. A file is a FASTA (with a `.fai` beside it, it is read
  from that), or a chromosome-length table with the tag `ft:cl`. Other tags are
  ignored. Names must not contain `#`.
- `--sr`: SyRI's `syri.out`, or plotsr's BEDPE (`achr astart aend bchr bstart
  bend type`), one per neighbouring pair in genomes-file order: the first
  compares genomes 1 and 2, the next 2 and 3.
- `--pair REF QRY FILE`: any other pair, named by genome.
- `--cfg`: plotsr's config file, read for `syncol`, `invcol`, `tracol` and
  `dupcol`. Unset, the colors are plotsr's defaults. Like plotsr, an inverted
  translocation or duplication takes the plain one's color.

SyRI pairs chromosomes by name, so homologous chromosomes are spelled alike in
every genome.

## Outputs

In the `-o` folder (default `syri_jbrowse`):

- `config.json`: one assembly per genome, the tracks below, and the stacked
  view it opens on
- `syri.paf`: every region as a PAF record between
  [PanSN](https://github.com/pangenome/PanSN-spec)-named sequences
  (`<genome>#1#<chromosome>`), its type in the tag `syri:Z:` and its color in
  `color:Z:`
- `syri_regions.bed`: the regions of each pair against the first genome, on
  that genome, a `query` column naming the other one
- `<name>.chrom.sizes`, one per genome

The config names every file relatively, so the folder moves as a whole.
