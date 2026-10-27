import os
import sys
import subprocess
import pysam
import argparse
import pandas as pd
import time
import array
import multiprocessing
from itertools import groupby
from sortedcontainers import SortedDict

def run_cmd(cmd):
    subprocess.run(cmd, shell=True, check=True)

def vcf_indexing(VCF, PATH):
    print('################################################' + '\n' + "Start VCF indexing!!!" + '\n' + '################################################')
    VCF_name = VCF.split('/')[-1].replace('.vcf', '')
    sorted_vcf = f'{PATH}/{VCF_name}.coordsorted.vcf'
    run_cmd(f"{{ grep '^#' {VCF} || true; grep -v '^#' {VCF} | sort -k1,1 -k2,2n; }} > {sorted_vcf}")
    run_cmd(f'bgzip -k -c -f {sorted_vcf} > {PATH}/{VCF_name}.vcf.gz')
    run_cmd(f'tabix -p vcf {PATH}/{VCF_name}.vcf.gz')
    os.remove(sorted_vcf)
    print('################################################' + '\n' + "Complete VCF indexing!!!" + '\n' + '################################################')

def make_bed(VCF, NAME, PATH):
    with open(f'{PATH}/{NAME}.vcf.bed', 'w') as note:
        with open(VCF, 'r') as vcf:
            for line in vcf:
                line = line.strip()
                if line.startswith('#'):
                    continue
                splitted = line.split('\t')
                Chr = str(splitted[0])
                Position_Start = str(int(splitted[1]) - 1)
                Position_End = str(int(splitted[1]))
                Ref = str(splitted[3])
                Alt = str(splitted[4])
                note.write(Chr + '\t' + Position_Start + '\t' + Position_End + '\t' + Ref + '\t' + Alt + '\n')
    print('################################################' + '\n' + "Complete VCF to BED!!!" + '\n' + '################################################')

def Split_BAM(BAM, PATH, NAME, THREAD):
    run_cmd(f'samtools view -@ {THREAD} -h -b {BAM} -L {PATH}/{NAME}.vcf.bed -U {PATH}/{NAME}.Normal.bam -o {PATH}/{NAME}.Mut.bam')

    pysam.sort("-@", f"{THREAD}", "-o", f'{PATH}/{NAME}.Mut.sorted.bam', f'{PATH}/{NAME}.Mut.bam')
    pysam.index("-@", f"{THREAD}", f'{PATH}/{NAME}.Mut.sorted.bam')
    os.remove(f'{PATH}/{NAME}.Mut.bam')

    print('################################################' + '\n' + "Complete splitting Mut & Normal BAM!!!" + '\n' + '################################################')

class CigarIntronMismatch(ValueError):
    pass


def Cigar_Check_NEW(CIGARLIST, N_lengths=None):
    Untied_list = [item for sublist in CIGARLIST
                   for item in (sublist if isinstance(sublist, list) else [sublist])]
    result = [(int(key), len(list(group))) for key, group in groupby(Untied_list)]

    if not N_lengths:
        if any(op == 3 for op, _ in result):
            raise CigarIntronMismatch(f"N operations present but no original intron lengths: {result}")
        return result

    fixed, ni = [], 0
    for op, n in result:
        if op == 3:
            for _ in range(n):
                if ni >= len(N_lengths):
                    raise CigarIntronMismatch(f"Too many N operations: {result} vs {N_lengths}")
                fixed.append((3, N_lengths[ni]))
                ni += 1
        else:
            fixed.append((op, n))
    if ni != len(N_lengths):
        raise CigarIntronMismatch(f"Too few N operations: {result} vs {N_lengths}")
    return fixed

def load_vcf(NAME, PATH):
    print("Loading VCF...")
    vcf = pd.read_csv(f"{PATH}/{NAME}.vcf.bed", sep="\t", header=None, low_memory=False)
    vcf.columns = ["Chromosome", "Start", "End", "Ref", "Alt"]

    vcf_dict = SortedDict()
    for _, row in vcf.iterrows():
        key = (str(row["Chromosome"]), row["Start"])
        if key not in vcf_dict:
            vcf_dict[key] = []
        vcf_dict[key].append((row["End"], row["Ref"], row["Alt"]))

    return vcf_dict

def query_vcf_range(vcf_data, chrom, bam_start, bam_end):
    if hasattr(vcf_data, 'irange'):
        return list(vcf_data.irange((chrom, bam_start), (chrom, bam_end)))
    return [k for k in vcf_data if k[0] == chrom and bam_start <= k[1] <= bam_end]

def parse_tag(tag_str):
    parts = tag_str.split(':', 2)
    if len(parts) != 3:
        return None
    name, type_code, value = parts
    try:
        if type_code == 'i':
            return (name, int(value))
        elif type_code == 'f':
            return (name, float(value))
        else:
            return (name, value)
    except (ValueError, TypeError):
        return (name, value)

def build_position_index(bam_read):
    CIGAR_string_Original = [i[0] for i in bam_read.cigartuples]
    CIGAR_number_Original = [i[1] for i in bam_read.cigartuples]
    CIGAR_list = sum(
        [list(str(CIGAR_string_Original[i]) * CIGAR_number_Original[i])
         if str(CIGAR_string_Original[i]) != '3'
         else [str(CIGAR_string_Original[i])]
         for i in range(0, len(CIGAR_string_Original))],
        []
    )

    N_pos = 0
    LIST_N_POS = []
    for i in CIGAR_string_Original:
        if i == 3:
            LIST_N_POS.append(N_pos)
            N_pos += 1
        else:
            N_pos += 1

    Original_Splicing_Cigar_number = []
    if len(LIST_N_POS) >= 1:
        Original_Splicing_Cigar_number = [CIGAR_number_Original[pos] for pos in LIST_N_POS]

    All_position_list = [bam_read.pos + 1]
    Last_position = All_position_list[-1]
    order = 0
    for i in bam_read.cigartuples:
        if i[0] == 1 and order == 0:
            Last_position = All_position_list[-1]
            All_position_list = [All_position_list[-1] - 1] * i[1]
        elif (i[0] == 0 and order == 0) or (i[0] == 2 and order == 0):
            All_position_list = list(range(Last_position, Last_position + i[1]))
            Last_position = All_position_list[-1] + 1
        elif i[0] == 4 and order == 0:
            All_position_list = list(range(Last_position - i[1], Last_position))
            Last_position = All_position_list[-1] + 1
        elif i[0] == 5 and order == 0:
            All_position_list = ['H'] * i[1]
        elif i[0] == 3 and order == 0:
            All_position_list = ['N'] + All_position_list
            Last_position = All_position_list[-1]
        elif i[0] == 1 and order != 0:
            All_position_list += [Last_position - 1] * i[1]
            Last_position = All_position_list[-1] + 1
        elif (i[0] == 0 and order != 0) or (i[0] == 2 and order != 0) or (i[0] == 4 and order != 0):
            All_position_list += list(range(Last_position, Last_position + i[1]))
            Last_position = All_position_list[-1] + 1
        elif i[0] == 5 and order != 0:
            All_position_list += ['H'] * i[1]
        elif i[0] == 3 and order != 0:
            Last_position = Last_position + i[1]
            All_position_list += ['N']
        order += 1

    Accumulated_position = 0
    for i in bam_read.cigartuples:
        if i[0] == 1:
            All_position_list[Accumulated_position: Accumulated_position + i[1]] = \
                [All_position_list[Accumulated_position: Accumulated_position + i[1]]]
            CIGAR_list[Accumulated_position: Accumulated_position + i[1]] = \
                [CIGAR_list[Accumulated_position: Accumulated_position + i[1]]]
            Accumulated_position += 1
        elif i[0] == 3:
            Accumulated_position += 1
        else:
            Accumulated_position += i[1]

    return CIGAR_list, All_position_list, Original_Splicing_Cigar_number

def compute_nm_md(ref_seq, query_seq, cigartuples):
    nm = 0
    md_parts = []
    match_run = 0
    ref_i = 0
    qry_i = 0

    for op, length in cigartuples:
        if op in (0, 7, 8):
            for _ in range(length):
                rbase = ref_seq[ref_i]
                qbase = query_seq[qry_i]
                if rbase.upper() == qbase.upper():
                    match_run += 1
                else:
                    nm += 1
                    md_parts.append(str(match_run))
                    md_parts.append(rbase.upper())
                    match_run = 0
                ref_i += 1
                qry_i += 1
        elif op == 1:
            nm += length
            qry_i += length
        elif op == 2:
            md_parts.append(str(match_run))
            md_parts.append('^' + ref_seq[ref_i:ref_i + length].upper())
            match_run = 0
            nm += length
            ref_i += length
        elif op == 3:
            ref_i += length
        elif op == 4:
            qry_i += length
        elif op == 5:
            pass

    md_parts.append(str(match_run))
    md_string = ''.join(md_parts)
    return nm, md_string


def update_nm_md_tag(bam_read, tags, ref_fasta):
    ref_seq = ref_fasta.fetch(bam_read.reference_name, bam_read.reference_start, bam_read.reference_end)
    nm, md = compute_nm_md(ref_seq, bam_read.query_sequence, bam_read.cigartuples)
    new_tags = []
    for name, value in tags:
        if name in ('NM', 'MD'):
            continue
        new_tags.append((name, value))
    new_tags.append(('NM', nm))
    new_tags.append(('MD', md))
    return new_tags

def Eraser(read_batch, output_file, THREAD, HEADER, vcf_data, ref_fasta_path):
    print('################################################' + '\n' + "Start Erasing step!!!" + '\n' + '################################################')
    HEADER = pysam.AlignmentHeader.from_dict(HEADER)
    output_bam = pysam.AlignmentFile(output_file, "wb", header=HEADER, threads=int(THREAD))
    ref_fasta = pysam.FastaFile(ref_fasta_path)
    mapping = {'1': 1, '2': -1, '5': -1}
    n_intron_skip = 0

    for read_str in read_batch:
        fields = read_str.split('\t')
        bam_read = pysam.AlignedSegment(header=HEADER)
        bam_read.query_name = fields[0]
        bam_read.query_sequence = fields[9]
        bam_read.flag = int(fields[1])
        bam_read.reference_name = fields[2]
        bam_read.reference_start = int(fields[3]) - 1
        bam_read.mapping_quality = int(fields[4])
        bam_read.cigarstring = fields[5]
        bam_read.next_reference_name = fields[6]
        bam_read.next_reference_start = int(fields[7]) - 1
        bam_read.template_length = int(fields[8])
        bam_read.query_qualities = pysam.qualitystring_to_array(fields[10])
        bam_read.tags = [t for t in (parse_tag(f) for f in fields[11:]) if t is not None]

        if bam_read.is_unmapped:
            output_bam.write(bam_read)
            continue

        if bam_read.cigartuples is None:
            output_bam.write(bam_read)
            continue

        bam_start = int(bam_read.reference_start)
        bam_end = int(bam_read.reference_end)
        chrom = bam_read.reference_name
        Pos_list, Ref, Alt = [], [], []

        for key in query_vcf_range(vcf_data, chrom, bam_start, bam_end):
            _, start = key
            for end, ref, alt in vcf_data[key]:
                if start <= bam_end and end >= bam_start:
                    Pos_list.append(start + 1)
                    Ref.append(ref)
                    Alt.append(alt)

        if len(Pos_list) == 0:
            output_bam.write(bam_read)
            continue

        CIGAR_list, All_position_list, Original_Splicing_Cigar_number = build_position_index(bam_read)
        cigar_dirty = False
        read_modified = False

        try:
            for idx in range(0, len(Pos_list)):
                if cigar_dirty:
                    CIGAR_list, All_position_list, Original_Splicing_Cigar_number = build_position_index(bam_read)
                    cigar_dirty = False

                NML = 0

                if len(Ref[idx]) == len(Alt[idx]):
                    if Pos_list[idx] not in All_position_list:
                        continue
                    TARGET = All_position_list.index(Pos_list[idx])
                    if CIGAR_list[TARGET] == '2':
                        continue
                    TARGET_list = CIGAR_list[0:TARGET]
                    NML = -len(list(filter(lambda x: isinstance(x, list), TARGET_list)))
                    TARGET_list = [item for sublist in TARGET_list for item in (sublist if isinstance(sublist, list) else [sublist])]
                    NML_list = list(filter(lambda x: (x != '0') and (x != '4'), TARGET_list))
                    if len(NML_list) > 0:
                        NML += sum(mapping.get(i, 0) for i in NML_list)
                    if '3' in TARGET_list:
                        NML += -TARGET_list.count('3')

                elif len(Ref[idx]) < len(Alt[idx]):
                    if [Pos_list[idx]] * (len(Alt[idx]) - 1) not in All_position_list or str(len(Alt[idx]) - len(Ref[idx])) + 'I' not in bam_read.cigarstring:
                        continue
                    TARGET = All_position_list.index(Pos_list[idx])
                    if CIGAR_list[TARGET] == '2':
                        continue
                    if set(CIGAR_list[TARGET + 1]) != {'1'}:
                        continue
                    TARGET_list = CIGAR_list[0:TARGET]
                    NML = -len(list(filter(lambda x: isinstance(x, list), TARGET_list)))
                    TARGET_list = [item for sublist in TARGET_list for item in (sublist if isinstance(sublist, list) else [sublist])]
                    NML_list = list(filter(lambda x: (x != '0') and (x != '4'), TARGET_list))
                    if len(NML_list) > 0:
                        NML += sum(mapping.get(i, 0) for i in NML_list)
                    if '3' in TARGET_list:
                        NML += -TARGET_list.count('3')

                elif len(Ref[idx]) > len(Alt[idx]):
                    if Pos_list[idx] not in All_position_list or str(len(Ref[idx]) - len(Alt[idx])) + 'D' not in bam_read.cigarstring:
                        continue
                    TARGET = All_position_list.index(Pos_list[idx]) + 1
                    if CIGAR_list[TARGET - 1] == '2':
                        continue
                    if set([item for sublist in CIGAR_list[TARGET: TARGET + len(Ref[idx]) - len(Alt[idx])] for item in (sublist if isinstance(sublist, list) else [sublist])]) != {'2'}:
                        continue
                    TARGET_list = CIGAR_list[0:TARGET]
                    NML = -len(list(filter(lambda x: isinstance(x, list), TARGET_list)))
                    TARGET_list = [item for sublist in TARGET_list for item in (sublist if isinstance(sublist, list) else [sublist])]
                    NML_list = list(filter(lambda x: (x != '0') and (x != '4'), TARGET_list))
                    if len(NML_list) > 0:
                        NML += sum(mapping.get(i, 0) for i in NML_list)
                    if '3' in TARGET_list:
                        NML += -TARGET_list.count('3')
                else:
                    continue

                Real_Index = NML

                if len(Ref[idx]) == len(Alt[idx]):
                    locus_idx = TARGET + Real_Index
                    Sequence = list(bam_read.query_sequence)
                    Qual = list(bam_read.query_qualities)
                    TPcode = bam_read.query_qualities.typecode
                    if ''.join(Sequence[locus_idx: locus_idx + len(Ref[idx])]) == Alt[idx]:
                        Sequence[locus_idx: locus_idx + len(Ref[idx])] = list(Ref[idx])
                        bam_read.query_sequence = ''.join(Sequence)
                        bam_read.query_qualities = array.array(TPcode, Qual)
                        read_modified = True

                elif len(Ref[idx]) < len(Alt[idx]) and len(CIGAR_list[TARGET + 1]) == len(Alt[idx]) - len(Ref[idx]):
                    locus_idx = TARGET + 1 + Real_Index
                    Sequence = list(bam_read.query_sequence)
                    Qual = list(bam_read.query_qualities)
                    TPcode = bam_read.query_qualities.typecode
                    ins_len = len(Alt[idx][len(Ref[idx]):])
                    del CIGAR_list[TARGET + 1]
                    new_cigartuples = Cigar_Check_NEW(CIGAR_list, Original_Splicing_Cigar_number)
                    del Sequence[locus_idx: locus_idx + ins_len]
                    bam_read.query_sequence = ''.join(Sequence)
                    bam_read.cigartuples = new_cigartuples
                    del Qual[locus_idx: locus_idx + ins_len]
                    bam_read.query_qualities = array.array(TPcode, Qual)
                    cigar_dirty = True
                    read_modified = True

                elif len(Ref[idx]) > len(Alt[idx]) and len(CIGAR_list[TARGET: TARGET + len(Ref[idx]) - len(Alt[idx])]) == len(Ref[idx]) - len(Alt[idx]):
                    locus_idx = TARGET + Real_Index
                    Sequence = list(bam_read.query_sequence)
                    Qual = list(bam_read.query_qualities)
                    TPcode = bam_read.query_qualities.typecode
                    del_len = len(Ref[idx]) - len(Alt[idx])
                    for j, base in enumerate(Ref[idx][len(Alt[idx]):]):
                        Sequence.insert(locus_idx + j, base)
                    avg_qual = round(sum(Qual) / len(Qual))
                    for j in range(del_len):
                        Qual.insert(locus_idx + j, avg_qual)
                    CIGAR_list[TARGET: TARGET + len(Ref[idx]) - len(Alt[idx])] = ['0'] * (len(Ref[idx]) - len(Alt[idx]))
                    new_cigartuples = Cigar_Check_NEW(CIGAR_list, Original_Splicing_Cigar_number)
                    bam_read.query_sequence = ''.join(Sequence)
                    bam_read.query_qualities = array.array(TPcode, Qual)
                    bam_read.cigartuples = new_cigartuples
                    cigar_dirty = True
                    read_modified = True

        except CigarIntronMismatch as e:
            n_intron_skip += 1
            if n_intron_skip <= 5:
                print(f"[warn] intron reassembly failed -> read left unmodified: "
                      f"{bam_read.query_name} {bam_read.cigarstring} ({e})", flush=True)

        tags = [t for t in (parse_tag(f) for f in fields[11:]) if t is not None]
        if read_modified:
            tags = update_nm_md_tag(bam_read, tags, ref_fasta)
        bam_read.tags = tags
        output_bam.write(bam_read)

    ref_fasta.close()
    output_bam.close()
    if n_intron_skip:
        print(f"[warn] reads passed unmodified after intron reassembly failure: {n_intron_skip} "
              f"({output_file}). Target variants may remain in these reads.", flush=True)

def main(MUT, PATH, NAME, THREAD, ref_fasta_path):
    print('################################################' + '\n' + "Start Main!!!" + '\n' + '################################################')
    Input_bam = pysam.AlignmentFile(MUT, "rb", threads=int(THREAD))
    Header_dict = Input_bam.header.to_dict()
    Read_list = [read.to_string() for read in Input_bam.fetch()]
    Input_bam.close()
    vcf_data = load_vcf(NAME, PATH)

    if THREAD > 1:
        chunk_size = len(Read_list) // THREAD

        if chunk_size >= THREAD:
            batches = [Read_list[chunk::THREAD] for chunk in range(THREAD)]
            vcf_data_plain = dict(vcf_data)

            with multiprocessing.Pool(THREAD) as pool:
                pool.starmap(Eraser, [(batch, f"{PATH}/Erased_tmp{num}.bam", THREAD, Header_dict, vcf_data_plain, ref_fasta_path) for num, batch in enumerate(batches)])

            run_cmd(f"samtools merge -@ {THREAD} -f -o {PATH}/Erased.bam {PATH}/Erased_tmp*")
            print('################################################ Complete Merging !!! ################################################', flush=True)
            run_cmd(f"rm -rf {PATH}/Erased_tmp*")
        else:
            print('################################################' + '\n' + "Start single processing!!!" + '\n' + '################################################')
            Eraser(Read_list, f"{PATH}/Erased.bam", THREAD, Header_dict, vcf_data, ref_fasta_path)
            print('################################################' + '\n' + "DONE !!!" + '\n' + '################################################')
    else:
        print('################################################' + '\n' + "Start single processing!!!" + '\n' + '################################################')
        Eraser(Read_list, f"{PATH}/Erased.bam", THREAD, Header_dict, vcf_data, ref_fasta_path)
        print('################################################' + '\n' + "DONE !!!" + '\n' + '################################################')

def headers_equal(bam_a, bam_b):
    IGNORE_KEYS = {'PG'}
    with pysam.AlignmentFile(bam_a, "rb") as a, pysam.AlignmentFile(bam_b, "rb") as b:
        da = {k: v for k, v in a.header.to_dict().items() if k not in IGNORE_KEYS}
        db = {k: v for k, v in b.header.to_dict().items() if k not in IGNORE_KEYS}
        return da == db

def is_declared_coordinate_sorted(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return bam.header.to_dict().get('HD', {}).get('SO') == 'coordinate'

def Merging(PATH, NAME, THREAD):
    print('################################################' + '\n' + 'Merging bam files', flush=True)
    run_cmd(f'samtools merge -@ {THREAD} -f -o {PATH}/{NAME}.Varaser.bam {PATH}/{NAME}.Normal.bam {PATH}/Erased.bam')
    print('DONE' + '\n' + '################################################', flush=True)

    if headers_equal(f'{PATH}/{NAME}.Normal.bam', f'{PATH}/{NAME}.Varaser.bam'):
        os.rename(f'{PATH}/{NAME}.Varaser.bam', f'{PATH}/{NAME}.Varaser.header.bam')
        print('################################################' + '\n' + "Header already matches source BAM - skipping reheader" + '\n' + '################################################', flush=True)
    else:
        run_cmd(f'samtools view -H {PATH}/{NAME}.Normal.bam > {PATH}/Header.txt')
        run_cmd(f'samtools reheader {PATH}/Header.txt {PATH}/{NAME}.Varaser.bam > {PATH}/{NAME}.Varaser.header.bam')

def sorting(PATH, NAME, THREAD, force_sort=False):
    src = f'{PATH}/{NAME}.Varaser.header.bam'
    dst = f'{PATH}/{NAME}.Varaser.sorted.bam'

    if not force_sort and is_declared_coordinate_sorted(src):
        os.rename(src, dst)
        print('################################################' + '\n' + "Input BAM already sorted - skipping full re-sort" + '\n' + '################################################', flush=True)
    else:
        run_cmd(f'samtools sort -@ {THREAD} {src} -o {dst}')

    run_cmd(f"samtools index {dst}")

def bamtofastq(TYPE, PATH, NAME, THREAD):
    if TYPE == 'bam':
        pass
    elif TYPE == 'fastq':
        run_cmd(f'samtools fastq -@ {THREAD} {PATH}/{NAME}.Varaser.sorted.bam -1 {PATH}/{NAME}_1.fastq.gz -2 {PATH}/{NAME}_2.fastq.gz')

def run():
    parser = argparse.ArgumentParser(description='Varaser V2')
    parser.add_argument("B", metavar='-B bam', help="input bam file")
    parser.add_argument("V", metavar='-V vcf', help="input vcf file")
    parser.add_argument("P", metavar='-P prefix', help="Set output prefix")
    parser.add_argument("L", metavar='-L prefix', help="Set output directory")

    parser.add_argument("-T", default='bam', type=str, dest='type', metavar='File type', help="Set your output file type (e.g. bam, fastq)")
    parser.add_argument("-N", default='1', type=str, dest='thread', metavar='thread', help="Set thread for run, default = 1")
    parser.add_argument("-R", "--reference", required=True,
                        dest='reference', metavar='reference fasta',
                        help="Reference FASTA used to recalculate NM/MD of modified reads (must match the BAM)")
    parser.add_argument("--force-sort", action='store_true', dest='force_sort',
                        help="Re-sort the final BAM even if it appears to be sorted already")
    args = parser.parse_args()

    input_bam = args.B
    Sample_name = args.P
    input_vcf = args.V
    WKDIR = args.L
    THREAD = int(args.thread)

    MutBAM = f'{WKDIR}/{Sample_name}.Mut.sorted.bam'
    os.makedirs(WKDIR, exist_ok=True)

    Start_time = time.time()

    vcf_indexing(VCF=input_vcf, PATH=WKDIR)
    make_bed(VCF=input_vcf, NAME=Sample_name, PATH=WKDIR)
    Split_BAM(BAM=input_bam, PATH=WKDIR, NAME=Sample_name, THREAD=THREAD)

    Erase_Start_time = time.time()
    main(MUT=MutBAM, PATH=WKDIR, NAME=Sample_name, THREAD=THREAD, ref_fasta_path=args.reference)
    Erase_End_time = time.time()

    Merging(PATH=WKDIR, NAME=Sample_name, THREAD=THREAD)
    sorting(PATH=WKDIR, NAME=Sample_name, THREAD=THREAD, force_sort=args.force_sort)

    End_time = time.time()
    Run_time = (End_time - Start_time) / 60
    Erasing_time = (Erase_End_time - Erase_Start_time) / 60

    with open(f'{WKDIR}/{Sample_name}.runtime.txt', 'w') as f:
        f.write(f'Run\t{Run_time:.3f}\n')
        f.write(f'Erase\t{Erasing_time:.3f}\n')

    run_cmd(
        f"rm -rf {WKDIR}/{Sample_name}.Mut.bam "
        f"{WKDIR}/{Sample_name}.vcf.bed "
        f"{WKDIR}/{Sample_name}.Normal.bam "
        f"{WKDIR}/{Sample_name}.Mut.sorted.bam "
        f"{WKDIR}/{Sample_name}.Mut.sorted.bam.bai "
        f"{WKDIR}/{Sample_name}.Varaser.header.bam "
        f"{WKDIR}/Erased.bam "
        f"{WKDIR}/Header.txt "
        f"{WKDIR}/{Sample_name}.Varaser.bam"
    )

if __name__ == '__main__':
    run()
