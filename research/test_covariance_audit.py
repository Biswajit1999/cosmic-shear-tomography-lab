import unittest
from covariance_audit import covariance_shrinkage_null, profile_constant_offset

class TestCovarianceAudit(unittest.TestCase):
    def test_diagonal(self):
        self.assertAlmostEqual(covariance_shrinkage_null([2.,3.],[[4.,0.],[0.,9.]])[0][1],2.)
    def test_correlated(self):
        out=covariance_shrinkage_null([1.,1.],[[1.,.5],[.5,1.]])
        self.assertAlmostEqual(out[0][1],4/3)
        self.assertAlmostEqual(out[-1][1],2.)
    def test_offset(self):
        r=profile_constant_offset([3.,3.,3.],[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        self.assertAlmostEqual(r["offset"],3.)
        self.assertAlmostEqual(r["sigma"],3**-.5)
        self.assertAlmostEqual(r["residual_chi2"],0.)
    def test_indefinite(self):
        with self.assertRaises(ValueError):
            profile_constant_offset([1.,2.],[[1.,2.],[2.,1.]])
if __name__=="__main__": unittest.main()
