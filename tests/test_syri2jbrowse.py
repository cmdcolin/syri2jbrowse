import json
import os
import tempfile
import unittest

from syri2jbrowse import main

# Two rows SyRI writes per region: the parent (column 10 '-') and an
# alignment under it, which is skipped. Coordinates are 1-based, inclusive;
# an inverted query runs end to start.
A_B = """\
Chr1\t1\t1000\t-\t-\tChr1\t1\t1000\tSYN1\t-\tSYN\t-
Chr1\t1\t1000\t-\t-\tChr1\t1\t1000\tSYNAL1\tSYN1\tSYNAL\t-
Chr1\t1001\t3000\t-\t-\tChr1\t3000\t1001\tINV1\t-\tINV\t-
Chr1\t3001\t3100\tA\tT\tChr1\t3001\t3001\tSNP1\t-\tSNP\t-
"""
# plotsr's BEDPE, which spells a translocation TRA
B_C = 'Chr1\t1\t1000\tChr2\t11\t1010\tTRA\n'
A_C = 'Chr1\t1\t4000\tChr1\t1\t4000\tSYN\n'


class Syri2JbrowseTest(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        d = self.dir.name
        files = {
            'A.chrlen': 'Chr1\t5000\n',
            'B.chrlen': 'Chr1\t5000\n',
            'C.chrlen': 'Chr1\t5000\nChr2\t2000\n',
            'genomes.txt': '#file\tname\ttags\nA.chrlen\tA\tft:cl\nB.chrlen\tB\tft:cl;lw:1.5\nC.chrlen\tC\tft:cl\n',
            'A_Bsyri.out': A_B,
            'B_C.bedpe': B_C,
            'A_C.bedpe': A_C,
            'base.cfg': 'syncol:#DEDEDE\ninvcol:#FFA500\n',
        }
        for name, text in files.items():
            with open(os.path.join(d, name), 'w') as fh:
                fh.write(text)
        self.out = os.path.join(d, 'out')

    def tearDown(self):
        self.dir.cleanup()

    def path(self, name):
        return os.path.join(self.dir.name, name)

    def run_tool(self, *extra):
        main(
            [
                '--genomes', self.path('genomes.txt'),
                '--sr', self.path('A_Bsyri.out'),
                '--sr', self.path('B_C.bedpe'),
                '-o', self.out,
                *extra,
            ]
        )
        with open(os.path.join(self.out, 'config.json')) as fh:
            return json.load(fh)

    def paf(self):
        with open(os.path.join(self.out, 'syri.paf')) as fh:
            return [line.rstrip('\n').split('\t') for line in fh]

    def test_regions_become_half_open_paf_records(self):
        self.run_tool()
        rows = self.paf()
        self.assertEqual(
            [r[:9] + r[12:13] for r in rows],
            [
                ['B#1#Chr1', '5000', '0', '1000', '+', 'A#1#Chr1', '5000', '0', '1000', 'syri:Z:SYN'],
                ['B#1#Chr1', '5000', '1000', '3000', '-', 'A#1#Chr1', '5000', '1000', '3000', 'syri:Z:INV'],
                ['C#1#Chr2', '2000', '10', '1010', '+', 'B#1#Chr1', '5000', '0', '1000', 'syri:Z:TRANS'],
            ],
        )
        self.assertEqual(rows[1][13], 'color:Z:#FFA500')

    def test_a_plotsr_config_recolors(self):
        self.run_tool('--cfg', self.path('base.cfg'))
        self.assertEqual(self.paf()[0][13], 'color:Z:#DEDEDE')

    def test_lanes_need_every_genome_against_the_first(self):
        config = self.run_tool()
        self.assertNotIn('displays', config['tracks'][0])
        config = self.run_tool('--pair', 'A', 'C', self.path('A_C.bedpe'))
        self.assertEqual(config['tracks'][0]['displays'][0]['rows'], {'sort': ['B', 'C']})

    def test_reference_rows_name_their_genome(self):
        config = self.run_tool('--pair', 'A', 'C', self.path('A_C.bedpe'))
        with open(os.path.join(self.out, 'syri_regions.bed')) as fh:
            rows = [line.split('\t') for line in fh.read().splitlines()[1:]]
        self.assertEqual([(r[1], r[3], r[9]) for r in rows], [('0', 'SYN', 'B'), ('0', 'SYN', 'C'), ('1000', 'INV', 'B')])
        color = config['tracks'][1]['displays'][0]['encoding']['color']
        self.assertEqual(color['scale']['domain'], ['SYN', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP'])
        self.assertEqual(color['scale']['range'][2], color['scale']['range'][3])

    def test_the_stack_hides_short_regions_as_plotsr_does(self):
        view = self.run_tool()['defaultSession']['views'][0]
        self.assertEqual(view['filter'], [{'field': 'length', 'gte': 10000}])
        view = self.run_tool('-s', '0')['defaultSession']['views'][0]
        self.assertNotIn('filter', view)

    def test_one_sr_per_neighbouring_pair(self):
        with self.assertRaises(SystemExit):
            main(['--genomes', self.path('genomes.txt'), '--sr', self.path('A_Bsyri.out'), '-o', self.out])


if __name__ == '__main__':
    unittest.main()
