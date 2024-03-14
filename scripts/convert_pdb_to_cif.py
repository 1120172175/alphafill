import os
from multiprocessing import Pool

import pandas as pd
from tqdm import tqdm
from Bio.PDB import MMCIFIO, PDBParser, PDBIO, MMCIFParser


SAVE_FOLDER = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/ECUST_coop/ESMFold_pred/cif'
INPUT_FOLDER = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/ECUST_coop/ESMFold_pred/pdb_processed'
TOOL_PATH = '/mnt/nas/ai-algorithm-data/liuyong/dataset/protein/PDB-REDO/tools/tools/pdb2cif'
CIF_TO_PDB_PATH = '/mnt/nas/ai-algorithm-data/liuyong/dataset/protein/PDB-REDO/tools/tools/cif2pdb'

def tranverse_folder(folder):
    filepath_list = []
    for root, _, files in os.walk(folder):
        for file in files:
            filepath_list.append(os.path.join(root, file))
    return filepath_list
    

def run_by_biopython(uniprot_id):
    if not isinstance(uniprot_id, str):
        return

    af_structure_path = f'{INPUT_FOLDER}/{uniprot_id}.pdb'
    if not os.path.exists(af_structure_path):
        print(f'AlphaFold structure not exists: {uniprot_id}')
        return
    
    save_path = f'{SAVE_FOLDER}/{uniprot_id}.cif'
    
    parser = PDBParser()
    structure = parser.get_structure("PDB_structure", af_structure_path)

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(save_path)


def run(uniprot_id):
    if not isinstance(uniprot_id, str):
        return

    af_structure_path = f'{INPUT_FOLDER}/{uniprot_id}.pdb'
    if not os.path.exists(af_structure_path):
        print(f'AlphaFold structure not exists: {uniprot_id}')
        return
    
    save_path = f'{SAVE_FOLDER}/{uniprot_id}.cif'
    
    command = f'{TOOL_PATH} {af_structure_path} {save_path}'
    os.system(command)




def main_multiprocess():
    list_path = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/rxneny_pro_all/uniprot_id_list.csv'
    uniprot_list = pd.read_csv(list_path)['uniprot_id'].values
    with Pool(30) as pool:
        for _ in tqdm(pool.imap(run, uniprot_list), total=len(uniprot_list)):
            pass


def main():
    subset_list_path = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/rxneny_pro_all/uniprot_id_list.csv'
    subset_uniprot_list = pd.read_csv(subset_list_path)['uniprot_id'].values

    for uniprot_id in tqdm(subset_uniprot_list):
        if not isinstance(uniprot_id, str):
            continue

        af_structure_path = f'/mnt/nas/ai-algorithm-data/liuyong/dataset/SynBio/enzyme-reaction-pairs/from_zt/3D-structures/alphafold_structures/{uniprot_id}.pdb'
        if not os.path.exists(af_structure_path):
            print(f'AlphaFold structure not exists: {uniprot_id}')
        save_path = f'/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/rxneny_pro_all/structures/{uniprot_id}.cif'

        running_command = f'{TOOL_PATH} {af_structure_path} {save_path}'
        os.system(running_command)


def main_ECUST_coop():
    pdb_path_list = tranverse_folder(INPUT_FOLDER)
    for pdb_path in tqdm(pdb_path_list):
        output_path = os.path.join(SAVE_FOLDER, os.path.basename(pdb_path).replace('.pdb', '.cif'))
        command = f'{TOOL_PATH} {pdb_path} {output_path}'
        os.system(command)


def tranverse_folder(folder):
    all_files = []
    for root, _, files in os.walk(folder):
        for file in files:
            all_files.append(os.path.join(root, file))
    return all_files


def main_cif_to_pdb():
    data_dir = '/gxr/liuyong/datasets/SynBio/af2/all_active_sites/cif'
    as_file_list = [each for each in tranverse_folder(data_dir) if '.cif' in each]
    as_file_list = [file for file in as_file_list if os.path.getsize(file) > 1000]
    print(f"Find {len(as_file_list)} structures.")

    for filepath in tqdm(as_file_list):
        parser = MMCIFParser()
        structure = parser.get_structure("protein", filepath)

        chain_names = [each.get_id() for each in structure.get_chains()]
        num_chains = len(chain_names)
        if num_chains > 1:
            print(os.path.basename(filepath).replace('.cif', ''), num_chains)
            print(chain_names)

        # filename = os.path.basename(filepath)
        # save_path = os.path.join(os.path.dirname(data_dir), 'pdb', filename.replace('.cif', '.pdb'))
        # io = PDBIO()
        # io.set_structure(structure)
        # io.save(save_path)


if __name__=="__main__":
    # main()
    # main_multiprocess()
    # main_ECUST_coop()
    main_cif_to_pdb()
