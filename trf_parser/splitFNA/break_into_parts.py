# !pip install biopython

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

for dolphin in ['candace', 'harish', 'ziheng']:
  records = list(SeqIO.parse("/home/ziheng/courses/CSE280/tool/splitFNA/data/compressedFNA_{}.fna".format(dolphin), "fasta"))
  part_records = []
  part_ix = 0
  curr_len = 0
  for record in tqdm(records):
    scaffold_info = re.search(r'scaffold(.*?),', record.description).group(1)
    scaffold_num, subscaffold_num = scaffold_info.split('_')[-2:]

    description = 'scaffold_{}_{}'.format(scaffold_num, subscaffold_num)
    new_record = record = SeqRecord(record.seq,
                                    description=description,
                                    id=record.id)
    part_records.append(new_record)

    curr_len += len(record.seq) + len(description) + len(record.id)
    if curr_len > 950000:
      with open('compressedFNA_{}_part_{}.fna'.format(dolphin, part_ix), 'w') as f:
        SeqIO.write(part_records, f, 'fasta')
      curr_len = 0
      part_ix += 1
      part_records = []

  with open('compressedFNA_{}_part_{}.fna'.format(dolphin, part_ix), 'w') as f:
      SeqIO.write(part_records, f, 'fasta')