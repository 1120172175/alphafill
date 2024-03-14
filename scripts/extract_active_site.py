import os
import string
from multiprocessing import Pool
from functools import partial
from collections import defaultdict

from Bio.PDB import MMCIFParser, PDBParser, PDBIO, NeighborSearch, MMCIFIO
import numpy as np
from tqdm import tqdm


def remove_far_residues(input_path, save_dir=None, cutoff=8.0, save_format='cif', preserve_ligand=True, remove_ion=False):
    try:
        if save_dir is None:
            save_dir = os.path.dirname(input_path)

        filename = os.path.basename(input_path)
        if '_transplant.cif' in filename:
            mark = '_transplant.cif'
        elif '_final.cif' in filename:
            mark = '_final.cif'
        save_path = os.path.join(save_dir, filename.replace(mark, f'.{save_format}'))

        if input_path.endswith('.cif'):
            parser = MMCIFParser()
            structure = parser.get_structure("protein", input_path)
        elif input_path.endswith('.pdb'):
            parser = PDBParser()
            structure = parser.get_structure("protein", input_path)

        structure = filter_ligands(structure, remove_ion=remove_ion)

        # 获取所有原子对象
        atoms = list(structure.get_atoms())

        ligand_atoms = []
        for residue in structure.get_residues():
            if residue.id[0] != ' ':
                ligand_atoms.extend(residue.get_atoms())

        # 创建NeighborSearch对象
        ns = NeighborSearch(atoms)

        # 找出与符合条件的小分子原子距离小于cutoff的氨基酸
        residues_to_keep = set()
        for ligand_atom in ligand_atoms:
            nearby_residues = ns.search(ligand_atom.coord, cutoff, level='R')
            for res in nearby_residues:
                residues_to_keep.add(res)

        # 删除距离大于cutoff的氨基酸
        for model in structure:
            for chain in model:
                for res in list(chain.get_residues()):
                    if res not in residues_to_keep:
                        chain.detach_child(res.id)

        if not preserve_ligand:
            structure = filter_ligands(structure, all_ligands=True)

        structure = reset_chain_id(structure)

        # 保存修改后的结构
        if save_format == 'cif':
            io = MMCIFIO()
        elif save_format == 'pdb':
            io = PDBIO()
        else:
            raise ValueError(f"Save format not supported: {save_format}")
        io.set_structure(structure)
        io.save(save_path)
    except Exception as e:
        print(e)
        return


def reset_chain_id(structure):
    for model in structure:
        for i, chain in enumerate(model):
            chain.id = string.ascii_uppercase[i]
    return structure


def filter_ligands(structure, all_ligands=False, remove_ion=False):
    if not all_ligands:
        ligand_residues_list = []
        for chain in structure.get_chains():
            residues = [residue for residue in chain]
            if len(residues) == 1:
                # ligand通常以单一的residue来表征
                ligand_residues_list.append(residues)
        
        ligand_cnt_dict = defaultdict(int)
        for residues in ligand_residues_list:
            resname = residues[0].resname
            ligand_cnt_dict[resname] += 1
        
        ligand_to_remove = set(['HOH', 'H2O'])
        for ligand_name, cnt in ligand_cnt_dict.items():
            if cnt >= 3:
                ligand_to_remove.add(ligand_name)
        
        chain_id_to_delete = []
        for model in structure:
            for chain in model:
                residues = [res for res in chain]
                if len(residues) > 1:
                    continue
                elif len(residues) == 1:
                    if residues[0].resname in ligand_to_remove:
                        chain_id_to_delete.append(chain.id)
                    elif remove_ion:
                        if len([_ for _ in residues[0].get_atoms()]) == 1:
                            chain_id_to_delete.append(chain.id)
                else:
                    raise ValueError(f"Error of len(residues): {len(residues)}")
    else:
        chain_id_to_delete = []
        for model in structure:
            for chain in model:
                residues = [res for res in chain]
                if len(residues) == 1:
                    chain_id_to_delete.append(chain.id)
    
    chain_id_to_delete = set(chain_id_to_delete)
    for model in structure:
        for chain_id in chain_id_to_delete:
            model.detach_child(chain_id)
        
    return structure
            


def main():
    # input_folder = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/output/rxneny_subset_5candidates/transplantion/'
    input_folder = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/output/rxneny_subset_all/transplantion'
    # input_folder = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/output/ECUST_coop/ESMFold_pred/transplantion_rerun'

    radius = 8
    # save_dir = '/gxr/liuyong/datasets/SynBio/af2/all_active_sites/pdb_10A'
    save_dir = '/mnt/nas/ai-algorithm-data/liuyong/dataset/SynBio/enzyme-reaction-pairs/from_zt/3D-structures/active_site_new_no_ion/pdb_8A'
    os.makedirs(save_dir, exist_ok=True)

    uniprot_id_list = []
    for subdir, dirs, files in tqdm(os.walk(input_folder)):
        for filename in files:
            if '_transplant.cif' in filename:
                uniprot_id_list.append(filename.replace('_transplant.cif', ''))
    
    input_list = []
    for uniprot_id in tqdm(uniprot_id_list):
        input_cif = f'{input_folder}/{uniprot_id}_transplant.cif'
        input_list.append(input_cif)
    
    running_func = partial(remove_far_residues, save_dir=save_dir, save_format='pdb', preserve_ligand=False, cutoff=radius, remove_ion=True)
    
    with Pool(40) as pool:
        for _ in tqdm(pool.imap(running_func, input_list), total=len(input_list)):
            pass



def main_pdb_redo():
    input_folder = '/mnt/nas/ai-algorithm-data/liuyong/dataset/protein/PDB-REDO/pdb-redo'
    pdb_id_list = []
    filepath_list = []
    print('Tranversing all files...')
    for subdir, dirs, files in tqdm(os.walk(input_folder)):
        for filename in files:
            if '_final.cif' in filename:
                pdb_id_list.append(filename.replace('_final.cif', ''))
                filepath_list.append(os.path.join(subdir, filename))
    
    print('Start extracting active site...')
    with Pool(30) as pool:
        for _ in tqdm(pool.imap(remove_far_residues, filepath_list), total=len(filepath_list)):
            pass
    

def move_result():
    input_folder = '/mnt/nas/ai-algorithm-data/liuyong/dataset/protein/PDB-REDO/pdb-redo'
    output_folder = '/mnt/nas/ai-algorithm-data/liuyong/dataset/protein/PDB-REDO/active_site'
    for subdir, dirs, files in tqdm(os.walk(input_folder)):
        for filename in files:
            if '_active_site.cif' in filename:
                filepath = os.path.join(subdir, filename)
                save_path = os.path.join(output_folder, filename)
                os.system(f'mv {filepath} {save_path}')


def remove_ligands():
    input_file = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/output/tmp/remove_ligands/B9MDV4_transplant.cif'
    save_dir = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/output/tmp/remove_ligands/'
    output_file = os.path.join(save_dir, os.path.basename(input_file).replace('_transplant.cif', '_trans_modified.pdb'))
    # print(output_file)
    remove_far_residues(input_file, save_dir=save_dir, save_format='pdb', preserve_ligand=False)

    # parser = MMCIFParser()
    # structure = parser.get_structure("protein", input_file)
    # structure = filter_ligands(structure)
    # structure = reset_chain_id(structure)
    # io = PDBIO()
    # io.set_structure(structure)
    # io.save(output_file)


def PaFS():
    trans_result = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/ECUST_coop/crystal/PaFS/trans_result/5er8_monomer_transplant.cif'
    save_path = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/ECUST_coop/crystal/PaFS/trans_result/5er8_monomer_trans_clean.pdb'

    parser = MMCIFParser()
    structure = parser.get_structure("protein", trans_result)
    structure = filter_ligands(structure)
    structure = reset_chain_id(structure)
    io = PDBIO()
    io.set_structure(structure)
    io.save(save_path)


def tranverse_folder(folder):
    all_files = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            all_files.append(os.path.join(root, file))
    return all_files


def rename():
    data_dir = '/gxr/liuyong/datasets/SynBio/af2/all_active_sites/pdb'
    file_list = tranverse_folder(data_dir)

    for filepath in tqdm(file_list):
        save_path = filepath.replace('_active_site.pdb', '.pdb')
        os.rename(filepath, save_path)


def delete_empty():
    # data_dir = '/gxr/liuyong/datasets/SynBio/af2/all_active_sites/pdb_10A'
    data_dir = '/mnt/nas/ai-algorithm-data/liuyong/dataset/SynBio/enzyme-reaction-pairs/from_zt/3D-structures/active_site_new_no_ion/pdb_8A'
    file_list = tranverse_folder(data_dir)
    for filepath in tqdm(file_list):
        if os.path.getsize(filepath) < 1000:
            os.remove(filepath)


if __name__ == "__main__":
    main()
    # main_pdb_redo()
    # move_result()
    # remove_ligands()
    # rename()
    delete_empty()
    # PaFS()