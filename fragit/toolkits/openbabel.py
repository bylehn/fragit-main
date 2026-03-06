"""
Copyright (C) 2011-2023 Casper Steinmann
"""
from __future__ import annotations
from typing import Iterable, Optional

from fragit.fragit_exceptions import OBNotFoundException

try:
    from openbabel import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")

from fragit.util import file_extension, Z2LABEL


class OBAtomWrapper:
    """Wraps openbabel.OBAtom, implementing AtomProtocol."""

    def __init__(self, ob_atom: openbabel.OBAtom):
        self._ob_atom = ob_atom

    def get_idx(self) -> int:
        return self._ob_atom.GetIdx()

    def get_id(self) -> int:
        return self._ob_atom.GetId()

    def get_atomic_num(self) -> int:
        return self._ob_atom.GetAtomicNum()

    def get_type(self) -> str:
        return self._ob_atom.GetType()

    def get_formal_charge(self) -> int:
        return self._ob_atom.GetFormalCharge()

    def set_formal_charge(self, charge: int) -> None:
        self._ob_atom.SetFormalCharge(charge)

    def get_implicit_h_count(self) -> int:
        return self._ob_atom.GetImplicitHCount()

    def set_implicit_h_count(self, count: int) -> None:
        self._ob_atom.SetImplicitHCount(count)

    def get_explicit_degree(self) -> int:
        try:
            return self._ob_atom.GetExplicitDegree()
        except AttributeError:
            return self._ob_atom.GetValence()

    def get_implicit_valence(self) -> int:
        try:
            return self._ob_atom.GetImplicitValence()
        except AttributeError:
            return self.get_explicit_degree() + self._ob_atom.GetImplicitHCount()

    def is_hbond_donor_h(self) -> bool:
        return self._ob_atom.IsHbondDonorH()

    def is_hbond_donor(self) -> bool:
        return self._ob_atom.IsHbondDonor()

    def is_hbond_acceptor(self) -> bool:
        return self._ob_atom.IsHbondAcceptor()

    def is_connected(self, other: OBAtomWrapper) -> bool:
        return self._ob_atom.IsConnected(other._ob_atom)

    def get_distance(self, other: OBAtomWrapper) -> float:
        return self._ob_atom.GetDistance(other._ob_atom)

    def get_angle(self, a1: OBAtomWrapper, a2: OBAtomWrapper) -> float:
        return self._ob_atom.GetAngle(a1._ob_atom, a2._ob_atom)

    def set_position(self, x: float, y: float, z: float) -> None:
        self._ob_atom.SetVector(x, y, z)

    def get_position(self) -> tuple[float, float, float]:
        return (self._ob_atom.GetX(), self._ob_atom.GetY(), self._ob_atom.GetZ())

    def iter_neighbors(self) -> Iterable[OBAtomWrapper]:
        for nb in openbabel.OBAtomAtomIter(self._ob_atom):
            yield OBAtomWrapper(nb)

    def __eq__(self, other) -> bool:
        if isinstance(other, OBAtomWrapper):
            return self._ob_atom == other._ob_atom
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self._ob_atom.GetIdx())


class OBBondWrapper:
    """Wraps openbabel.OBBond, implementing BondProtocol."""

    def __init__(self, ob_bond: openbabel.OBBond):
        self._ob_bond = ob_bond

    def get_bond_order(self) -> int:
        return self._ob_bond.GetBondOrder()


class OBResidueWrapper:
    """Wraps an OpenBabel residue, implementing ResidueProtocol."""

    def __init__(self, ob_residue):
        self._ob_residue = ob_residue

    def get_name(self) -> str:
        return self._ob_residue.GetName()

    def get_atom_id(self, atom: OBAtomWrapper) -> str:
        return self._ob_residue.GetAtomID(atom._ob_atom)

    def get_num_atoms(self) -> int:
        return self._ob_residue.GetNumAtoms()

    def iter_atoms(self) -> Iterable[OBAtomWrapper]:
        for atom in openbabel.OBResidueAtomIter(self._ob_residue):
            yield OBAtomWrapper(atom)


class OBMolWrapper:
    """Wraps openbabel.OBMol, implementing MoleculeProtocol."""

    def __init__(self, ob_mol: openbabel.OBMol):
        self._ob_mol = ob_mol

    def num_atoms(self) -> int:
        return self._ob_mol.NumAtoms()

    def get_atom(self, idx: int) -> OBAtomWrapper:
        return OBAtomWrapper(self._ob_mol.GetAtom(idx))

    def delete_atom(self, atom: OBAtomWrapper) -> None:
        self._ob_mol.DeleteAtom(atom._ob_atom)

    def add_atom(self, atom: OBAtomWrapper) -> bool:
        return self._ob_mol.AddAtom(atom._ob_atom)

    def duplicate_atom(self, source: OBAtomWrapper) -> OBAtomWrapper:
        new_ob_atom = openbabel.OBAtom()
        new_ob_atom.Duplicate(source._ob_atom)
        return OBAtomWrapper(new_ob_atom)

    def get_bond(self, idx1: int, idx2: int) -> Optional[OBBondWrapper]:
        bond = self._ob_mol.GetBond(idx1, idx2)
        if bond is None:
            return None
        return OBBondWrapper(bond)

    def delete_bond(self, bond: OBBondWrapper) -> None:
        self._ob_mol.DeleteBond(bond._ob_bond)

    def find_children(self, a1: int, a2: int) -> list[int]:
        tmp = openbabel.vectorInt()
        self._ob_mol.FindChildren(tmp, a2, a1)
        return [value for value in tmp]

    def add_hydrogens_to_atom(self, atom: OBAtomWrapper) -> bool:
        return self._ob_mol.AddHydrogens(atom._ob_atom)

    def iter_residues(self) -> Iterable[OBResidueWrapper]:
        for residue in openbabel.OBResidueIter(self._ob_mol):
            yield OBResidueWrapper(residue)


class OBSmartsMatcher:
    """SMARTS pattern matcher using openbabel.OBSmartsPattern."""

    def __init__(self):
        self._pattern = openbabel.OBSmartsPattern()

    def match(self, mol: OBMolWrapper, pattern: str) -> list[tuple[int, ...]]:
        self._pattern.Init(pattern)
        self._pattern.Match(mol._ob_mol)
        return [m for m in self._pattern.GetUMapList()]


class OBChargeModelWrapper:
    """Charge model wrapper using openbabel.OBChargeModel."""

    def __init__(self, model_name: str):
        self._charge_model = openbabel.OBChargeModel.FindType(model_name)
        if self._charge_model is None:
            raise ValueError("Charge model '{}' not found.".format(model_name))

    def compute_charges(self, mol: OBMolWrapper) -> bool:
        return self._charge_model.ComputeCharges(mol._ob_mol)

    def get_partial_charges(self) -> list[float]:
        return list(self._charge_model.GetPartialCharges())


class OBToolkitAdapter:
    """OpenBabel implementation of ToolkitAdapterProtocol."""

    def read_molecule(self, filename: str) -> OBMolWrapper:
        file_format = file_extension(filename)[1:]
        ob_mol = openbabel.OBMol()
        conversion = openbabel.OBConversion()
        conversion.SetInFormat(file_format)
        conversion.ReadFile(ob_mol, filename)
        return OBMolWrapper(ob_mol)

    def make_smarts_matcher(self) -> OBSmartsMatcher:
        return OBSmartsMatcher()

    def make_charge_model(self, model_name: str) -> OBChargeModelWrapper:
        return OBChargeModelWrapper(model_name)

    def make_empty_atom(self) -> OBAtomWrapper:
        return OBAtomWrapper(openbabel.OBAtom())


class Molecule(object):
    """Legacy molecule wrapper kept for backwards compatibility."""

    def __init__(self, filename):
        self._setup(filename)

    def _setup(self, filename):
        self._adapter = OBToolkitAdapter()
        self._mol_wrapper = self._adapter.read_molecule(filename)
        self._compute_partial_atomic_charges()
        self._setup_smarts_pattern()

    def _compute_partial_atomic_charges(self):
        self._charge_model = self._adapter.make_charge_model("mmff94")
        self._charge_model.compute_charges(self._mol_wrapper)

    def _setup_smarts_pattern(self):
        self._matcher = self._adapter.make_smarts_matcher()

    def is_ok(self):
        return self.get_atom_count() > 1

    def match_pattern(self, pattern_to_match):
        return self._matcher.match(self._mol_wrapper, pattern_to_match)

    def get_partial_atom_charges(self):
        return self._charge_model.get_partial_charges()

    def get_total_charge(self) -> int:
        return int(sum(self.get_partial_atom_charges()))

    @staticmethod
    def get_element_symbol(atom_index: int) -> str:
        return Z2LABEL[atom_index]

    def get_atom_count(self):
        return self._mol_wrapper.num_atoms()
