"""SyRI's structural annotations as a JBrowse 2 config.

Takes what plotsr takes, a genomes file and one annotation per neighbouring
pair, and writes beside a config.json:

- syri.paf, every region as a PAF record between PanSN-named sequences, its
  type in the tag syri:Z: and that type's color in color:Z:
- syri_regions.bed, the regions of every pair whose reference is the first
  genome, placed on that genome, a `query` column naming the other one
- <name>.chrom.sizes per genome
"""

import argparse
import gzip
import json
import os
import sys

__version__ = '0.1.0'

TYPES = ['SYN', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP']

# plotsr's config keys and its default colors (config/base.cfg). plotsr draws
# an inverted translocation or duplication as the plain one, and so does this.
PLOTSR_COLORS = {
    'syncol': '#CCCCCC',
    'invcol': '#FFA500',
    'tracol': '#9ACD32',
    'dupcol': '#00BBFF',
}
COLOR_KEY = {
    'SYN': 'syncol',
    'INV': 'invcol',
    'TRANS': 'tracol',
    'INVTR': 'tracol',
    'DUP': 'dupcol',
    'INVDP': 'dupcol',
}
# plotsr's BEDPE input names a translocation TRA
TYPE_ALIASES = {'TRA': 'TRANS'}


def read_cfg(path):
    colors = dict(PLOTSR_COLORS)
    if path:
        with open(path) as fh:
            for line in fh:
                key, _, value = line.strip().partition(':')
                if key in colors and value:
                    colors[key] = value.strip()
    return colors


def parse_tags(text):
    return dict(t.split(':', 1) for t in text.split(';') if ':' in t)


def read_genomes(path):
    base = os.path.dirname(os.path.abspath(path))
    genomes = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2:
                sys.exit(f'{path}: each genome line is <file>\\t<name>[\\t<tags>]: {line.strip()}')
            file, name = fields[0], fields[1]
            tags = parse_tags(fields[2]) if len(fields) > 2 else {}
            resolved = file if os.path.exists(file) else os.path.join(base, file)
            genomes.append({'name': name, 'file': resolved, 'ft': tags.get('ft', 'fa')})
    names = [g['name'] for g in genomes]
    if len(set(names)) != len(names):
        sys.exit(f'{path}: genome names repeat')
    bad = [n for n in names if '#' in n]
    if bad:
        sys.exit(f'{path}: a genome name cannot hold "#", which PanSN sequence names split on: {bad[0]}')
    return genomes


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def fasta_lengths(path):
    lengths = {}
    name = None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def chrom_lengths(genome):
    path = genome['file']
    if genome['ft'] == 'cl':
        table = path
    elif os.path.exists(path + '.fai'):
        table = path + '.fai'
    else:
        return fasta_lengths(path)
    with open(table) as fh:
        return {f[0]: int(f[1]) for f in (line.split('\t') for line in fh if line.strip())}


def read_regions(path):
    """The top-level structural regions of a syri.out, or of plotsr's BEDPE
    (achr astart aend bchr bstart bend type), 1-based and inclusive as SyRI
    writes them, as half-open intervals."""
    found = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) >= 11:
                kind, parent = f[10], f[9]
                if parent != '-':
                    continue
                ref = (f[0], f[1], f[2])
                qry = (f[5], f[6], f[7])
            elif len(f) >= 7:
                kind = f[6]
                ref = (f[0], f[1], f[2])
                qry = (f[3], f[4], f[5])
            else:
                sys.exit(f'{path}: neither a syri.out row nor a BEDPE row: {line.strip()}')
            kind = TYPE_ALIASES.get(kind, kind)
            if kind not in COLOR_KEY:
                continue
            rstart, rend = sorted((int(ref[1]), int(ref[2])))
            qstart, qend = sorted((int(qry[1]), int(qry[2])))
            found.append(
                {
                    'type': kind,
                    'ref': (ref[0], rstart - 1, rend),
                    'qry': (qry[0], qstart - 1, qend),
                    'strand': '-' if kind.startswith('INV') else '+',
                }
            )
    return found


def paf_lines(found, reference, query, lengths, colors):
    for r in found:
        rchrom, rstart, rend = r['ref']
        qchrom, qstart, qend = r['qry']
        for genome, chrom in ((reference, rchrom), (query, qchrom)):
            if chrom not in lengths[genome]:
                sys.exit(f'{chrom} is in the SyRI output for {reference} and {query} but not in {genome}\'s sequences')
        spans = (rend - rstart, qend - qstart)
        yield (
            f'{query}#1#{qchrom}\t{lengths[query][qchrom]}\t{qstart}\t{qend}\t{r["strand"]}'
            f'\t{reference}#1#{rchrom}\t{lengths[reference][rchrom]}\t{rstart}\t{rend}'
            f'\t{min(spans)}\t{max(spans)}\t60'
            f'\tsyri:Z:{r["type"]}\tcolor:Z:{colors[COLOR_KEY[r["type"]]]}\n'
        )


BED_COLUMNS = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'thickStart',
    'thickEnd',
    'itemRgb',
    'query',
]


def rgb(hex_color):
    h = hex_color.lstrip('#')
    return ','.join(str(int(h[i : i + 2], 16)) for i in (0, 2, 4))


def bed_rows(found, query, colors):
    for r in found:
        chrom, start, end = r['ref']
        yield (chrom, start, end, r['type'], r['strand'], rgb(colors[COLOR_KEY[r['type']]]), query)


def color_scale(colors):
    return {'domain': TYPES, 'range': [colors[COLOR_KEY[t]].lower() for t in TYPES]}


def config(names, colors, lane_pairs, has_regions, min_size):
    reference = names[0]
    syri_color = {'field': 'syri', 'scale': {'domain': TYPES}}
    synteny_track = {
        'type': 'SyntenyTrack',
        'trackId': 'syri',
        'name': 'SyRI regions',
        'assemblyNames': names,
        'adapter': {
            'type': 'MultiGenomePAFAdapter',
            'uri': 'syri.paf',
            'attributeColumns': ['syri', 'color'],
        },
    }
    if lane_pairs:
        synteny_track['displays'] = [
            {
                'type': 'MultiWaySyntenyDisplay',
                'rows': {'sort': names[1:]},
                'ribbons': {'encoding': {'color': syri_color, 'opacity': 0.7}},
            }
        ]
    tracks = [synteny_track]
    if has_regions:
        tracks.append(
            {
                'type': 'FeatureTrack',
                'trackId': f'syri_regions_on_{reference}',
                'name': f'SyRI regions on {reference}, by genome',
                'assemblyNames': [reference],
                'adapter': {
                    'type': 'BedAdapter',
                    'uri': 'syri_regions.bed',
                    'columnNames': BED_COLUMNS,
                    'disableGeneHeuristic': True,
                },
                'displays': [
                    {
                        'type': 'LinearMultiRowFeatureDisplay',
                        'rows': {'field': 'query', 'sort': names[1:]},
                        'encoding': {
                            'color': {
                                'field': 'name',
                                'scale': color_scale(colors),
                                'legend': {'title': 'syri'},
                            }
                        },
                    }
                ],
            }
        )
    view = {
        'type': 'LinearSyntenyView',
        'views': [{'assembly': n} for n in names],
        'tracks': ['syri'],
        'encoding': {'color': syri_color, 'opacity': 0.9},
        **({'filter': [{'field': 'length', 'gte': min_size}]} if min_size > 0 else {}),
        'mark': {'interpolate': 'bump'},
        'render': {'antialias': False},
        'collapseEmptyRows': True,
        'levelHeights': [150] * (len(names) - 1),
    }
    return {
        '$schema': 'https://jbrowse.org/jb2/schema/v5/config.json',
        'assemblies': [
            {
                'name': n,
                'sequence': {
                    'type': 'ReferenceSequenceTrack',
                    'trackId': f'{n}-ReferenceSequenceTrack',
                    'adapter': {'type': 'ChromSizesAdapter', 'uri': f'{n}.chrom.sizes'},
                },
            }
            for n in names
        ],
        'tracks': tracks,
        'defaultSession': {'name': 'SyRI regions', 'views': [view]},
    }


def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog='syri2jbrowse',
        description=(
            "Write a JBrowse 2 config showing SyRI's syntenic, inverted, translocated and "
            'duplicated regions between genomes. Takes the inputs plotsr takes.'
        ),
    )
    parser.add_argument(
        '--sr',
        action='append',
        default=[],
        metavar='FILE',
        help='syri.out (or plotsr BEDPE) of one neighbouring pair: the first --sr compares genomes 1 and 2, the next 2 and 3, as plotsr reads them',
    )
    parser.add_argument(
        '--pair',
        action='append',
        default=[],
        nargs=3,
        metavar=('REF', 'QRY', 'FILE'),
        help="any other pair, named by genome: give one for each genome against the first, and JBrowse draws every genome as lanes in the first one's coordinates",
    )
    parser.add_argument('--genomes', required=True, help="plotsr's genomes file: <file>\\t<name>[\\t<tags>], ft:cl marking a chromosome-length file")
    parser.add_argument('--cfg', help="plotsr's config file, read for syncol, invcol, tracol and dupcol")
    parser.add_argument('-s', '--min-size', type=int, default=10000, help='the stacked view hides regions shorter than this, as plotsr -s does (default 10000; 0 shows all)')
    parser.add_argument('-o', '--out', default='syri_jbrowse', help='output directory (default syri_jbrowse)')
    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    genomes = read_genomes(args.genomes)
    names = [g['name'] for g in genomes]
    if len(args.sr) != len(names) - 1 and args.sr:
        sys.exit(f'{len(names)} genomes take {len(names) - 1} --sr files, one per neighbouring pair; got {len(args.sr)}')
    pairs = [(names[i], names[i + 1], path) for i, path in enumerate(args.sr)]
    for ref, qry, path in args.pair:
        for n in (ref, qry):
            if n not in names:
                sys.exit(f'--pair names {n}, which {args.genomes} does not list')
        pairs.append((ref, qry, path))
    if not pairs:
        sys.exit('give at least one --sr or --pair')
    seen = set()
    for ref, qry, _ in pairs:
        if (ref, qry) in seen:
            sys.exit(f'{ref} against {qry} is given twice')
        seen.add((ref, qry))

    colors = read_cfg(args.cfg)
    lengths = {g['name']: chrom_lengths(g) for g in genomes}
    os.makedirs(args.out, exist_ok=True)
    for name, table in lengths.items():
        with open(os.path.join(args.out, f'{name}.chrom.sizes'), 'w') as fh:
            fh.writelines(f'{chrom}\t{size}\n' for chrom, size in table.items())

    reference = names[0]
    bed = []
    with open(os.path.join(args.out, 'syri.paf'), 'w') as paf:
        for ref, qry, path in pairs:
            found = read_regions(path)
            paf.writelines(paf_lines(found, ref, qry, lengths, colors))
            if ref == reference:
                bed.extend(bed_rows(found, qry, colors))
            print(f'{ref} vs {qry}: {len(found)} regions', file=sys.stderr)
    if bed:
        bed.sort(key=lambda r: (r[0], r[1], r[2]))
        with open(os.path.join(args.out, 'syri_regions.bed'), 'w') as fh:
            fh.write('#' + '\t'.join(BED_COLUMNS) + '\n')
            for chrom, start, end, kind, strand, item_rgb, qry in bed:
                fh.write(f'{chrom}\t{start}\t{end}\t{kind}\t0\t{strand}\t{start}\t{end}\t{item_rgb}\t{qry}\n')

    against_reference = {qry for ref, qry, _ in pairs if ref == reference}
    lanes = against_reference >= set(names[1:])
    with open(os.path.join(args.out, 'config.json'), 'w') as fh:
        json.dump(config(names, colors, lanes, bool(bed), args.min_size), fh, indent=2)
        fh.write('\n')
    print(f'wrote {os.path.join(args.out, "config.json")}', file=sys.stderr)
