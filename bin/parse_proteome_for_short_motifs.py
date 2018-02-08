import sys
import os
import re

def parse_fasta_file(input_file):
    seq_uid_dict = {}
    with open(input_file, 'r') as in_fh:
        for entry in in_fh.read().split('\n>'):
            protein_id = entry.strip().split('\n')[0].replace('>', '')
            protein_seq_list = ''.join(entry.split('\n')[1:])
            seq_uid_dict[protein_id] = protein_seq_list
    return seq_uid_dict

def get_doisorder_coordinates(disorder_uid_dict):
    disorder_coordinate_dict = {}
    for entry in disorder_uid_dict.keys():
        uid = entry.split('|')[1]
        start = entry.split('|')[2].split('_')[-2]
        end = entry.split('|')[2].split('_')[-1]
        disorder_coordinate_dict.setdefault(
            uid, []).append("%s\t%s\t%s\t%s\tDisorder\t" % (
            uid, disorder_uid_dict[entry], start, end))
    return disorder_coordinate_dict

def parse_complete_proteome(proteome_data_dict, disorder_coordinate_dict):
    order_coordinate_dict = {}
    for prot_id in proteome_data_dict.keys():
        if prot_id in disorder_coordinate_dict.keys():
            prot_seq = proteome_data_dict[prot_id]
            prot_len = len(prot_seq)
            for i, each_disorder in enumerate(disorder_coordinate_dict[prot_id]):
                start_coordinate = int(each_disorder.split('\t')[2])
                end_coordinate = int(each_disorder.split('\t')[3])
                if prot_id not in order_coordinate_dict.keys():
                    order_seq = prot_seq[1:start_coordinate]
                    order_coordinate_dict.setdefault(
                    prot_id, []).append("%s\t%s\t%s\t%s\tOrder\t" % (
                    prot_id, order_seq, 1, start_coordinate))
                elif i == len(disorder_coordinate_dict[prot_id]):
                    last_coordinate = int(disorder_coordinate_dict[
                        prot_id][i-1].split('\t')[3])
                    order_seq = prot_seq[last_coordinate:prot_len]
                    order_coordinate_dict.setdefault(
                    prot_id, []).append("%s\t%s\t%s\t%s\tOrder\t" % (
                    prot_id, order_seq, last_coordinate, prot_len))
                else:
                    last_coordinate = int(disorder_coordinate_dict[
                        prot_id][i-1].split('\t')[3])
                    order_seq = prot_seq[last_coordinate:start_coordinate]
                    order_coordinate_dict.setdefault(
                    prot_id, []).append("%s\t%s\t%s\t%s\tOrder\t" % (
                    prot_id, order_seq, last_coordinate, start_coordinate))
    return order_coordinate_dict
                    
def create_out_file(out_file, disorder_coordinate_dict,
                    order_coordinate_dict, motif_list):
    with open(out_file, 'w') as out_fh:
        header = '\t'.join(['UID', 'FastaSeq', 'Start',
                            'End', 'Order-Disorder', 'SQ', 'TQ', 'QS', 'QT'])
        out_fh.write(header+'\n')
        for each_id in disorder_coordinate_dict.keys():
            for order_entry in order_coordinate_dict[each_id]:
                order_seq = order_entry.split('\t')[1]
                motif_order_list = check_motif_occurance(order_seq, motif_list)
                out_fh.write(order_entry+'\t'.join(map(str, motif_order_list))+'\n')
            for disorder_entry in disorder_coordinate_dict[each_id]:
                disorder_seq = disorder_entry.split('\t')[1]
                motif_disorder_list = check_motif_occurance(disorder_seq, motif_list)
                out_fh.write(disorder_entry+'\t'.join(map(str, motif_order_list))+'\n')

def check_motif_occurance(sequence_entry, motif_list):
    motif_occur_list = []
    for motif in motif_list:
        if motif in sequence_entry:
            motif_occur_list.append(len(list(
                enumerate(re.finditer(motif, sequence_entry)))))
        else:
            motif_occur_list.append(0)
    return motif_occur_list

if __name__ == '__main__':
    proteome_file = sys.argv[1]
    disorder_file = sys.argv[2]
    out_path = sys.argv[3]
    motif_list = list(sys.argv[4].split(',')) #['SQ', 'TQ', 'QS', 'TS']
    disorder_data_dict = parse_fasta_file(disorder_file)
    disorder_coordinate_dict = get_doisorder_coordinates(disorder_data_dict)
    order_data_dict = parse_fasta_file(proteome_file)
    order_coordinate_dict = parse_complete_proteome(order_data_dict, disorder_coordinate_dict)
    create_out_file(out_path, disorder_coordinate_dict, order_coordinate_dict, motif_list)
     