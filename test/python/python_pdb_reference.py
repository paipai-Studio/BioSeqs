#!/usr/bin/env python3
"""Python reference script for Bio.PDB module comparison."""

from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue, Atom
import numpy as np

def create_test_structure():
    """Create a test structure for comparison."""
    structure = Structure.Structure("test")
    model = Model.Model(1)
    structure.add(model)
    
    chain_a = Chain.Chain('A')
    model.add(chain_a)
    
    # Add residues with atoms
    for i, resname in enumerate(['ALA', 'SER', 'LYS']):
        residue = Residue.Residue((' ', i+1, ' '), resname, 'A')
        
        # Add some atoms
        if resname == 'ALA':
            residue.add(Atom.Atom('N', [0.0, 0.0, 0.0], 1.0, 0.0, ' ', 'N', 1, 'N'))
            residue.add(Atom.Atom('CA', [1.5, 0.0, 0.0], 1.0, 0.0, ' ', 'C', 2, 'C'))
            residue.add(Atom.Atom('C', [3.0, 0.0, 0.0], 1.0, 0.0, ' ', 'C', 3, 'C'))
        elif resname == 'SER':
            residue.add(Atom.Atom('N', [4.5, 0.0, 0.0], 1.0, 0.0, ' ', 'N', 4, 'N'))
            residue.add(Atom.Atom('CA', [6.0, 0.0, 0.0], 1.0, 0.0, ' ', 'C', 5, 'C'))
            residue.add(Atom.Atom('CB', [6.0, 1.5, 0.0], 1.0, 0.0, ' ', 'C', 6, 'C'))
            residue.add(Atom.Atom('OG', [6.0, 3.0, 0.0], 1.0, 0.0, ' ', 'O', 7, 'O'))
        elif resname == 'LYS':
            residue.add(Atom.Atom('N', [7.5, 0.0, 0.0], 1.0, 0.0, ' ', 'N', 8, 'N'))
            residue.add(Atom.Atom('CA', [9.0, 0.0, 0.0], 1.0, 0.0, ' ', 'C', 9, 'C'))
        
        chain_a.add(residue)
    
    return structure

def test_structure_stats(structure):
    """Test structure statistics."""
    print("=== Structure Statistics ===")
    print(f"Number of models: {len(structure)}")
    
    for model in structure:
        print(f"  Model {model.id}")
        print(f"    Number of chains: {len(model)}")
        
        for chain in model:
            print(f"      Chain {chain.id}")
            print(f"        Number of residues: {len(chain)}")
            
            for residue in chain:
                print(f"          Residue {residue.resname} {residue.get_id()[1]}")
                print(f"            Number of atoms: {len(residue)}")
                
                for atom in residue:
                    print(f"              Atom {atom.name}: ({atom.coord[0]:.3f}, {atom.coord[1]:.3f}, {atom.coord[2]:.3f})")

def test_atom_distance(structure):
    """Test atom distance calculation."""
    print("\n=== Atom Distance Calculation ===")
    
    # Get CA atoms from first two residues
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.name == 'CA':
                        atoms.append(atom)
    
    if len(atoms) >= 2:
        dist = np.linalg.norm(np.array(atoms[0].coord) - np.array(atoms[1].coord))
        print(f"Distance between CA atoms: {dist:.4f}")
        
        # Expected: 4.5 Angstroms (6.0 - 1.5)
        print(f"Expected: 4.5000")

def test_center_of_mass(structure):
    """Test center of mass calculation."""
    print("\n=== Center of Mass Calculation ===")
    
    for model in structure:
        for chain in model:
            for residue in chain:
                coords = np.array([atom.coord for atom in residue])
                com = coords.mean(axis=0)
                print(f"Residue {residue.resname} COM: ({com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f})")

def test_pdb_io():
    """Test PDB file I/O."""
    print("\n=== PDB I/O Test ===")
    
    pdb_content = """ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1      3.000   0.000   0.000  1.00  0.00           C
ATOM      4  N   SER A   2      4.500   0.000   0.000  1.00  0.00           N
ATOM      5  CA  SER A   2      6.000   0.000   0.000  1.00  0.00           C
END
"""
    
    # Write to temp file
    with open("/tmp/test.pdb", "w") as f:
        f.write(pdb_content)
    
    # Parse with Bio.PDB
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("test", "/tmp/test.pdb")
    
    print(f"Parsed structure has {len(structure)} models")
    
    atom_count = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                atom_count += len(residue)
    
    print(f"Total atoms: {atom_count}")
    
    # Write back
    io = PDBIO()
    io.set_structure(structure)
    io.save("/tmp/test_out.pdb")
    
    with open("/tmp/test_out.pdb") as f:
        written_content = f.read()
    
    print("\nWritten PDB content:")
    print(written_content)

def test_structure_traversal(structure):
    """Test structure traversal."""
    print("\n=== Structure Traversal ===")
    
    # Depth-first traversal
    all_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    all_atoms.append(atom)
    
    print(f"Total atoms via traversal: {len(all_atoms)}")
    
    # Get CA atoms only
    ca_atoms = [atom for atom in all_atoms if atom.name == 'CA']
    print(f"CA atoms: {len(ca_atoms)}")
    for ca in ca_atoms:
        print(f"  CA in {ca.parent.resname} {ca.parent.get_id()[1]}")

if __name__ == "__main__":
    structure = create_test_structure()
    
    test_structure_stats(structure)
    test_atom_distance(structure)
    test_center_of_mass(structure)
    test_structure_traversal(structure)
    test_pdb_io()
    
    print("\n=== Reference Output Complete ===")