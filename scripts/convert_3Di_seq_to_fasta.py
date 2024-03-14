

def convert(input_db, save_path):
    raw_data = open(input_db, 'r').readlines()
    raw_data = [line.split('\t')[:3] for line in raw_data]

    content_to_save = []
    for each in raw_data:
        desc = each[0]
        uniprot_id = desc.split('.')[0]
        struc_seq = each[2]
        content_to_save.append(f'>{uniprot_id}\n')
        content_to_save.append(f'{struc_seq}\n')
    
    with open(save_path, 'w') as f:
        f.writelines(content_to_save)

    print(f'Saved to {save_path}')


if __name__ == "__main__":
    # input_db = '/mnt/nas/ai-algorithm-data/liuyong/dataset/SynBio/enzyme-reaction-pairs/from_zt/3D-structures/3Di_seq_db/db'
    # save_path = '/mnt/nas/ai-algorithm-data/liuyong/dataset/SynBio/enzyme-reaction-pairs/from_zt/3D-structures/3Di_seq_db/db.fasta'

    input_db = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/output/rxneny_subset_5candidates/active_site_3Di_seq/cutoff_10/db'
    save_path = '/mnt/nas/ai-algorithm-data/liuyong/code/SynBio/alphafill/output/rxneny_subset_5candidates/active_site_3Di_seq/cutoff_10/db.fasta'
    convert(input_db, save_path)
