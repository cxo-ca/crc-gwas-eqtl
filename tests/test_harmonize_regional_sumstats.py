import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from harmonize_regional_sumstats import align, deduplicate, key, normalize, reasons


def row(**changes):
    result = dict(locus='L1', chr='1', pos='100', effect_allele='A', other_allele='C',
                  beta='0.2', se='0.05', p='0.8', build='GRCh38', hm_code='10',
                  source_rsid='rs1', source_row=1, gene_id='ENSG1.2', source_gene_id='ENSG1',
                  tissue='Colon', sdY='0.98', source_variant='chr1_100_C_A')
    result.update(changes)
    return result


class HarmonizationTests(unittest.TestCase):
    def test_swap_changes_beta_only(self):
        g, e = row(), row(effect_allele='C', other_allele='A', beta='-0.3')
        z = align(g, e)
        self.assertEqual(z['allele_flipped'], 'TRUE')
        self.assertEqual(float(z['beta_eqtl_aligned']), 0.3)
        self.assertEqual((z['se_eqtl'], z['p_eqtl']), ('0.05', '0.8'))
        self.assertEqual(e['beta'], '-0.3')
        self.assertEqual(align(g, row())['allele_flipped'], 'FALSE')

    def test_no_rsid_only_or_complement_matching(self):
        for e in [row(pos='101'), row(effect_allele='T', other_allele='G'), row(build='GRCh37')]:
            with self.assertRaises(ValueError):
                align(row(), e)
        self.assertNotEqual(key(row()), key(row(other_allele='G')))

    def test_invalid_and_ambiguous_reason_overlap(self):
        why = reasons(row(effect_allele='G', other_allele='C', se='-nan'), 'eqtl')
        self.assertEqual(set(why), {'invalid_se', 'palindromic_snp'})
        for se in ['0', '-1', 'inf', '1e-300', '1e300']:
            self.assertIn('invalid_se', reasons(row(se=se), 'eqtl'))
        self.assertIn('invalid_p', reasons(row(p='0'), 'gwas'))
        self.assertIn('non_snv', reasons(row(effect_allele='AT'), 'eqtl'))
        self.assertIn('unharmonized_gwas', reasons(row(hm_code='15'), 'gwas'))
        self.assertEqual(reasons(row(p='1'), 'gwas'), [])  # no significance filtering

    def test_alias_duplicates_preserve_ids(self):
        kept, removed = deduplicate([row(), row(source_rsid='rs2', source_row=2)], 'eqtl')
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0]['source_rsid'], 'rs1;rs2')
        self.assertEqual(removed[0][1:], (['duplicate_rsid_alias'], 1))

    def test_conflicting_duplicate_quarantines_all(self):
        for change in [dict(beta='0.3'), dict(se='0.1'), dict(p='0.1')]:
            kept, removed = deduplicate([row(), row(source_row=2, **change)], 'eqtl')
            self.assertEqual(kept, [])
            self.assertEqual(len(removed), 2)
            self.assertTrue(all(r[1] == ['conflicting_duplicate'] for r in removed))

    def test_gene_and_tissue_are_distinct_groups(self):
        kept, removed = deduplicate([row(), row(gene_id='ENSG2.1'), row(tissue='Other')], 'eqtl')
        self.assertEqual((len(kept), len(removed)), (3, 0))

    def test_build_region_variant_and_sd_guards(self):
        regions = {'L1': dict(chr='1', start='50', end='150', build='GRCh38')}
        sdy = {('ENSG1.2', 'Colon'): 0.98}
        self.assertEqual(normalize(row(chr='chr1'), 'eqtl', regions, sdy)['chr'], '1')
        for change in [dict(build='GRCh37'), dict(pos='151'), dict(chr='2'),
                       dict(source_variant='chr1_100_A_C'), dict(sdY='1'), dict(source_gene_id='ENSG2')]:
            with self.assertRaises(ValueError):
                normalize(row(**change), 'eqtl', regions, sdy)


if __name__ == '__main__':
    unittest.main()
