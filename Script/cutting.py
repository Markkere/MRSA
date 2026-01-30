from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def cutting(file,contig,start,end,desired_file_name):
    with open(file) as org_file:
        id_list = []
        extracted_seq = []
        start_pos = int(start) -1 
        end_pos = int(end) + 1
        for record in SeqIO.parse(org_file,'fasta'):
            name_id = record.id
            if contig in name_id:
                desire_seq = record.seq[start_pos:end_pos]
                new_name_id = f"{contig}_{start}_{end}"
                id_list.append(new_name_id)
                extracted_seq.append(desire_seq)
            else:
                id_list.append(new_name_id)
                extracted_seq.append('Not found')
    final_seq = SeqRecord(extracted_seq[0],id = id_list[0],description='')
    with open(desired_file_name,'w') as des_name:
        SeqIO.write(final_seq,des_name,'fasta')    
        
cutting('/Users/marker/Desktop/MRSA_handover/Genome/SGH_N121_BC02_genomic.fna','contig_1',1949564,1990621,'BC02_orfX_sccmec_check.fasta')





    