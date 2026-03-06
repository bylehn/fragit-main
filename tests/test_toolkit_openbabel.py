"""
Copyright (C) 2026 Casper Steinmann
"""
from __future__ import annotations

import unittest

from fragit.toolkits.openbabel import OBToolkitAdapter
from fragit.toolkits import get_toolkit

WATER_CLUSTER = "tests/watercluster4.xyz"   # 4 H2O = 12 atoms


class TestOBToolkitAdapter(unittest.TestCase):

    def setUp(self):
        self.adapter = OBToolkitAdapter()
        self.mol = self.adapter.read_molecule(WATER_CLUSTER)

    def test_read_molecule_atom_count(self):
        """read_molecule returns a MoleculeProtocol with the correct number of atoms."""
        self.assertEqual(self.mol.num_atoms(), 12)

    def test_atom_properties(self):
        """OBAtomWrapper exposes protocol-compatible properties with sensible values."""
        atom = self.mol.get_atom(1)   # 1-based, first atom (an oxygen)
        self.assertIsInstance(atom.get_atomic_num(), int)
        self.assertGreater(atom.get_atomic_num(), 0)
        x, y, z = atom.get_position()
        self.assertIsInstance(x, float)
        self.assertIsInstance(y, float)
        self.assertIsInstance(z, float)

    def test_smarts_matcher_oxygen(self):
        """OBSmartsMatcher finds the correct number of oxygen atoms via SMARTS."""
        matcher = self.adapter.make_smarts_matcher()
        # watercluster4 has 4 water molecules → 4 oxygens
        matches = matcher.match(self.mol, "[#8]")
        self.assertEqual(len(matches), 4)

    def test_charge_model_partial_charges_length(self):
        """OBChargeModelWrapper.get_partial_charges returns one value per atom."""
        charge_model = self.adapter.make_charge_model("mmff94")
        charge_model.compute_charges(self.mol)
        charges = charge_model.get_partial_charges()
        self.assertEqual(len(charges), self.mol.num_atoms())

    def test_make_empty_atom_has_zero_atomic_num(self):
        """make_empty_atom returns an AtomProtocol with atomic number 0 (unset)."""
        atom = self.adapter.make_empty_atom()
        self.assertEqual(atom.get_atomic_num(), 0)


class TestGetToolkitFactory(unittest.TestCase):

    def test_get_toolkit_returns_adapter(self):
        """get_toolkit() returns an object satisfying ToolkitAdapterProtocol."""
        toolkit = get_toolkit("openbabel")
        self.assertTrue(hasattr(toolkit, "read_molecule"))
        self.assertTrue(hasattr(toolkit, "make_smarts_matcher"))
        self.assertTrue(hasattr(toolkit, "make_charge_model"))
        self.assertTrue(hasattr(toolkit, "make_empty_atom"))

    def test_get_toolkit_unknown_name_raises(self):
        """get_toolkit() raises ImportError for an unrecognised toolkit name."""
        with self.assertRaises(ImportError):
            get_toolkit("nonexistent_toolkit")


def suite():
    s = unittest.TestSuite()
    s.addTest(unittest.makeSuite(TestOBToolkitAdapter))
    s.addTest(unittest.makeSuite(TestGetToolkitFactory))
    return s


if __name__ == "__main__":
    unittest.main()
