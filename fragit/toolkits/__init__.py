"""
Copyright (C) 2026 Casper Steinmann
"""
import os
from typing import Optional

from fragit.toolkits.base import (
    AtomProtocol,
    BondProtocol,
    ResidueProtocol,
    MoleculeProtocol,
    SmartsMatcherProtocol,
    ChargeModelProtocol,
    ToolkitAdapterProtocol,
)

_toolkit_instance = None


def get_toolkit(name: Optional[str] = None) -> ToolkitAdapterProtocol:
    global _toolkit_instance
    if _toolkit_instance is not None and name is None:
        return _toolkit_instance

    toolkit_name = name or os.environ.get("FRAGIT_TOOLKIT", "openbabel")

    if toolkit_name == "openbabel":
        from fragit.toolkits.openbabel import OBToolkitAdapter
        _toolkit_instance = OBToolkitAdapter()
        return _toolkit_instance
    else:
        raise ImportError("Unknown toolkit: {!r}. Available: 'openbabel'".format(toolkit_name))
