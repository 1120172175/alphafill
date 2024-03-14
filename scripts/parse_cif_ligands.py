import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem

def cif_to_smiles(cif_path):
    doc = gemmi.cif.read_file(cif_path)

    smiles_list = []

    # 遍历CIF文件中的所有数据块
    for block in doc:
        chemcomp = gemmi.make_chemcomp_from_block(block)
        # print(chemcomp)
        # mol = chemcomp.make_molecule()
        # 查看mol这个instance的所有属性
        print(dir(chemcomp))
        break
        
        # 遍历所有化学组分
        for chem_comp in block.find_loop("_ligand.id"):
            print('sss')
            comp_id = chem_comp.str(0)
            struct = gemmi.make_chemcomp_structure(comp_id, block)

            # 转换为RDKit分子
            mol = struct.make_molecule()
            rdkit_mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))

            # 添加氢并计算2D坐标
            rdkit_mol = Chem.AddHs(rdkit_mol)
            AllChem.Compute2DCoords(rdkit_mol)

            # 生成SMILES字符串
            smiles = Chem.MolToSmiles(rdkit_mol)
            smiles_list.append(smiles)

    return smiles_list

# 示例使用
cif_path = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/af-ligands.cif'  # 替换为您的CIF文件路径
smiles_list = cif_to_smiles(cif_path)
print(smiles_list)
for smiles in smiles_list:
    print(smiles)
    break

