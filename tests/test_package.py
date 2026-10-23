import unittest

import cichlid_data
import cichlid_orthofinder


class TestPackage(unittest.TestCase):
    def test_coordinator_in_pinned_range(self):
        major, minor = (int(x) for x in cichlid_data.__version__.split(".")[:2])
        self.assertEqual((major, minor), (0, 14))

    def test_version(self):
        self.assertTrue(cichlid_orthofinder.__version__)


if __name__ == "__main__":
    unittest.main()
