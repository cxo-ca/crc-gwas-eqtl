"""Offline checks of the independent verifier, using saved synthetic R output."""
import math
import sys
import unittest
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import report_real_coloc as report


class ColocPosteriorTests(unittest.TestCase):
    def test_saved_synthetic_r_posteriors(self):
        # Frozen synthetic output; its sessionInfo records R 4.6.1 / coloc 5.2.3.
        # Read the saved fits rather than regenerating expectations with Python.
        demo = ROOT / 'results/demo'
        blocks = defaultdict(list)
        for row in report.rows(demo / 'harmonized.tsv'):
            blocks[report.group(row)].append(row)
        expected_groups = {
            ('SIMULATED_LOCUS', gene, 'SIMULATED_COLON')
            for gene in ('SIM_SHARED', 'SIM_DISTINCT')
        }
        self.assertEqual(set(blocks), expected_groups)
        stats = {}
        for key, block in blocks.items():
            self.assertEqual(len(block), 240)
            self.assertEqual(len({row['key'] for row in block}), 240)
            gwas = [report.labf(float(r['beta_gwas']), float(r['se_gwas']), 0.2)
                    for r in block]
            eqtl = [report.labf(float(r['beta_eqtl_aligned']), float(r['se_eqtl']),
                                0.15 * float(r['sdY'])) for r in block]
            stats[key] = (report.logsum(gwas), report.logsum(eqtl),
                          report.logsum([a + b for a, b in zip(gwas, eqtl)]))
        fits = list(report.rows(demo / 'coloc_sensitivity.tsv'))
        self.assertEqual(len(fits), 6)
        self.assertEqual({(report.group(r), float(r['p12'])) for r in fits},
                         {(key, prior) for key in expected_groups
                          for prior in (1e-6, 1e-5, 1e-4)})
        for fit in fits:
            key = report.group(fit)
            self.assertEqual(int(fit['nsnps']), len(blocks[key]))
            self.assertEqual((float(fit['p1']), float(fit['p2'])), (1e-4, 1e-4))
            actual = report.independent_posteriors(
                stats[key], *(float(fit[p]) for p in ('p1', 'p2', 'p12')))
            self.assertEqual(len(actual), 5)
            self.assertAlmostEqual(math.fsum(actual), 1.0, delta=1e-12)
            for hypothesis, value in enumerate(actual):
                with self.subTest(group=key, p12=fit['p12'], hypothesis=hypothesis):
                    expected = float(fit[f'PP_H{hypothesis}'])
                    self.assertTrue(math.isfinite(value) and 0 <= value <= 1)
                    # Relative tolerance also protects the tiny H0-H2 values.
                    self.assertTrue(math.isclose(value, expected, rel_tol=1e-10, abs_tol=0),
                                    (value, expected))

    def test_direct_hypothesis_enumeration(self):
        # Two SNPs with ABF1=(2,3), ABF2=(5,7): H3 has only 2*7 + 3*5.
        # One SNP cannot support distinct causal variants, so its H3 is zero.
        cases = [((math.log(5), math.log(12), math.log(31)),
                  [1, 0.05, 0.24, 0.0058, 0.031]),
                 ((math.log(2), math.log(5), math.log(10)),
                  [1, 0.02, 0.10, 0, 0.01])]
        for stats, weights in cases:
            with self.subTest(stats=stats):
                actual = report.independent_posteriors(stats, 0.01, 0.02, 0.001)
                self.assertEqual(len(actual), len(weights))
                for value, weight in zip(actual, weights):
                    self.assertAlmostEqual(value, weight / math.fsum(weights), delta=1e-14)


if __name__ == '__main__':
    unittest.main()
