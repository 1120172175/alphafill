import os
from multiprocessing import Pool

from tqdm import tqdm
import pandas as pd


def read_txt(path):
    with open(path, "r") as f:
        data = [each.strip() for each in f.readlines()]
    return data


def tranverse_folder(folder):
    filepath_list = []
    for root, _, files in os.walk(folder):
        for file in files:
            filepath_list.append(os.path.join(root, file))
    return filepath_list


def main(job_id):
    input_folder = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/rxneny_pro_subset/structures'
    uniprot_id_list = []
    for subdir, dirs, files in os.walk(input_folder):
        for filename in files:
            if '.cif' in filename:
                uniprot_id_list.append(filename.split('.')[0])
    
    chunk_size = len(uniprot_id_list) // 3
    if job_id == 0:
        input_list = uniprot_id_list[:chunk_size]
    elif job_id == 1:
        input_list = uniprot_id_list[chunk_size:chunk_size*2]
    elif job_id == 2:
        input_list = uniprot_id_list[chunk_size*2:]

    
    output_folder = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/output/rxneny_subset'
    for uniprot_id in tqdm(input_list):
        input_structure_path = f'{input_folder}/{uniprot_id}.cif'
        output_path = f'{output_folder}/{uniprot_id}_transplant.cif'
        command = f'alphafill process {input_structure_path} {output_path} --pdb-fasta=/mnt/nas/ai-algorithm-data/liuyong/dataset/protein/AlphaFill/pdbredo_seqdb.txt --pdb-dir=/mnt/nas/ai-algorithm-data/liuyong/dataset/protein/PDB-REDO/pdb-redo --blast-report-limit=5'
        os.system(command) 


def main_ECUST(rerun=False):
    input_folder = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/ECUST_coop/ESMFold_pred/cif'
    input_list = tranverse_folder(input_folder)

    if rerun:
        rerun_list = read_txt('/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/ECUST_coop/ESMFold_pred/rerun_proteins.txt')
        rerun_list = [each.replace('_active_site.cif', '.cif') for each in rerun_list]
        input_list = [os.path.join(input_folder, filename) for filename in rerun_list]
    
    output_folder = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/output/ECUST_coop/ESMFold_pred/transplantion'
    if rerun:
        output_folder = output_folder + '_rerun'
    os.makedirs(output_folder, exist_ok=True)

    blast_limit = 20
    if rerun:
        blast_limit = 200

    command_list = []
    for input_structure_path in tqdm(input_list):
        output_path = f'{output_folder}/{os.path.basename(input_structure_path).replace(".cif", "_transplant.cif")}'
        command = f'alphafill process {input_structure_path} {output_path} --pdb-fasta=/mnt/nas/ai-algorithm-data/liuyong/dataset/protein/AlphaFill/pdbredo_seqdb.txt --pdb-dir=/mnt/nas/ai-algorithm-data/liuyong/dataset/protein/PDB-REDO/pdb-redo --blast-report-limit={blast_limit}'
        command_list.append(command)
    
    with Pool(40) as p:
        for _ in tqdm(p.imap(os.system, command_list), total=len(command_list)):
            pass


def main_multiprocess(job_id=-1):
    input_folder = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/rxneny_pro_all/structures'
    uniprot_id_list = []
    for subdir, dirs, files in os.walk(input_folder):
        for filename in files:
            if '.cif' in filename:
                uniprot_id_list.append(filename.split('.')[0])
    # uniprot_id_list = pd.read_csv('/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/input/rxneny_pro_subset/rerun_uniprot_id_list.csv')['uniport_id'].tolist()

    chunk_size = len(uniprot_id_list) // 2
    if job_id == 0:
        uniprot_id_list = uniprot_id_list[:chunk_size]
    elif job_id == 1:
        uniprot_id_list = uniprot_id_list[chunk_size:]

    command_list = []
    output_folder = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/output/rxneny_subset_all/transplantion/'
    for uniprot_id in tqdm(uniprot_id_list):
        input_structure_path = f'{input_folder}/{uniprot_id}.cif'
        output_path = f'{output_folder}/{uniprot_id}_transplant.cif'
        command = f'alphafill process {input_structure_path} {output_path} --pdb-fasta=/mnt/nas/ai-algorithm-data/liuyong/dataset/protein/AlphaFill/pdbredo_seqdb.txt --pdb-dir=/mnt/nas/ai-algorithm-data/liuyong/dataset/protein/PDB-REDO/pdb-redo --blast-report-limit=5'
        command_list.append(command)
    
    with Pool(40) as p:
        for _ in tqdm(p.imap(os.system, command_list), total=len(command_list)):
            pass


if __name__=='__main__':
    # job_id = 2
    # main(job_id)

    # main_multiprocess(1)
    main_ECUST(rerun=True)
