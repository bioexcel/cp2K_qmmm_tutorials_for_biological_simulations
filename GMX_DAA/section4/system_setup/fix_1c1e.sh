# Split water molecules, protein and ligand ENH
grep "HOH" 1c1e.pdb > 1c1e.xray_waters.pdb
grep -v "HOH" 1c1e.pdb | grep -v "ENH" > 1c1e.protein.pdb
grep "ENH" 1c1e.pdb > 1c1e.ligand.pdb
