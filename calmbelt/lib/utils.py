from Bio.Blast import NCBIXML
import pandas as pd
import numpy as np
import os

# -------------- for alignment ----------------------------------------------------------

def hsp_to_snps(hsp, person):
    q = hsp.query  # reference
    ms = hsp.match
    s = hsp.sbjct  # patient
    start = hsp.query_start
    end = hsp.query_end

    snps_list = []

    for i, m in enumerate(ms):
        if m != '|':
            pos = i + start
            #             snps_list += ['{}{}{}'.format(q[i],pos,s[i])]
            snps_list += [[person, q[i], pos, s[i], "{}_{}".format(start, end), end - start + 1]]
    return snps_list


def start_to_end(hsp, person):
    start = hsp.query_start
    end = hsp.query_end
    return [person, "{}_{}".format(start, end), end - start + 1]


def start_end_list_func(url):
    result = []
    start_end = []
    snps_list = []
    start_end_list = []
    match = []

    for record in NCBIXML.parse(open(url)):
        print(record)
        if record.alignments:
            for i, align in enumerate(record.alignments):
                person = str(align).split('|')[1].split('_')[2]
                match += [person]

                for hsp in align.hsps:
                    snps_list += hsp_to_snps(hsp, person)
                    start_end_list += [start_to_end(hsp, person)]

                result += snps_list
                start_end += start_end_list
                snps_list = []
                start_end_list = []

    return result, start_end


# ไม่อ่านชื่่อ sample ในไฟล์
def start_end_list_func_built_tree(url):
    result = []
    start_end = []
    snps_list = []
    start_end_list = []
    match = []

    for record in NCBIXML.parse(open(url)):
        print(record)
        if record.alignments:
            for i, align in enumerate(record.alignments):
                person = str(align)
                match += [person]

                for hsp in align.hsps:
                    snps_list += hsp_to_snps(hsp, person)
                    start_end_list += [start_to_end(hsp, person)]

                result += snps_list
                start_end += start_end_list
                snps_list = []
                start_end_list = []

    return result, start_end


# -------------- for delete insertion ---------------------------------------------------

def drop_query(snps_max_df, sample_id, start_end):
    snps_id = snps_max_df[snps_max_df['id'] == sample_id]
    indexNames = snps_max_df[(snps_max_df['start_end'] == start_end) & (snps_max_df['query'] == '-')].index
    snps_max_df.drop(indexNames, inplace=True)
    return snps_max_df


def del_insertion(id_del, snps_df):
    # ลบ - และเลื่อนตำแหน่งสำหรับ snps ในตำแหน่งถัดๆไป

    for sample_id in id_del:
        start_end_del = set(snps_df[(snps_df['id'] == sample_id) & (snps_df['query'] == '-')]['start_end'])
        for sed in start_end_del:
            new_pos = []
            count = 1
            i = 0

            pos_del = sorted(list(
                snps_df[(snps_df['id'] == sample_id) & (snps_df['query'] == '-') & (snps_df['start_end'] == sed)][
                    'position']))
            snp_del = snps_df[(snps_df['id'] == sample_id) & (snps_df['start_end'] == sed)].sort_values(
                ['position']).copy()
            pos_all = sorted(list(snp_del['position']))

            for a in pos_all:
                old_pos = snp_del[snp_del['position'] == a]['position'].values[0]
                if a <= pos_del[i]:  # ไม่เปลี่ยน
                    new_pos += [old_pos]
                else:
                    if pos_del[i] == pos_del[-1]:  # a>pos_del[i] และ pos_del[i]เป็นตัวสุดท้าย
                        new_pos += [old_pos - count]
                    else:
                        if a < pos_del[i + 1]:  # a>pos_del[i] และ pos_del[i]ไมใช่ตัวสุดท้าย
                            new_pos += [old_pos - count]
                        else:
                            count += 1
                            i += 1
                            new_pos += [old_pos - count]

            snps_df.loc[(snps_df['id'] == sample_id) & (snps_df['start_end'] == sed), 'position'] = new_pos
            snps_df = drop_query(snps_df, sample_id, sed)
    return snps_df


# -------------- สำหรับหาตำแหน่่งทั้งหมดใน genome แล้วใส่ 1 -------------------------------------

def drop_data(snps_max_df, sample_id, intersec_set, start_end):
    snps_id = snps_max_df[snps_max_df['id'] == sample_id]
    position = list(snps_id[snps_id['start_end'] == start_end]['position'])
    for pos in position:
        if pos in intersec_set:
            indexNames = snps_max_df[(snps_max_df['start_end'] == start_end) & (snps_max_df['position'] == pos)].index
            snps_max_df.drop(indexNames, inplace=True)
    return snps_max_df


def length_max_snps(align_array, id_list, startstop_df, snps_max_df):
    # ตัด snps ตำแหน่งที่ซ้ำกันออก (โดยดูสายทั้งหมดจาก startstop_df )
    sum_length = []
    for i, sample_id in enumerate(id_list):  # id_list
        start_end_set = set()
        se_df = startstop_df[startstop_df.index == sample_id].sort_values('length', ascending=False)
        se = pd.unique(se_df['start_end'])  # list(set(se_df['start_end']))
        for j in se:
            (a, b) = tuple(map(int, j.split('_')))
            #         print(a,b)
            if (start_end_set) & set(range(a, b + 1)) == set():  # สายไม่ซ้ำกัน
                align_array[i][a - 1:b] = 1
                start_end_set = start_end_set.union(set(range(a, b + 1)))
            else:  # สายซ้ำกัน
                intersec_set = (start_end_set) & set(range(a, b + 1))
                new_pos = set(range(a, b + 1)) - intersec_set
                if new_pos != set():
                    for k in new_pos:
                        align_array[i][k - 1] = 1

                snps_max_df = drop_data(snps_max_df, sample_id, intersec_set, j)
                start_end_set = start_end_set.union(set(range(a, b + 1)))

        sum_length += [[sample_id, len(start_end_set)]]
    return sum_length, snps_max_df, align_array


def length_max_snps2(align_array, id_list, startstop_df, snps_max_df):  # สำหรับหา position insertion
    # ตัด snps ตำแหน่งที่ซ้ำกันออก (โดยดูสายทั้งหมดจาก startstop_df )
    sum_length = []
    choose_pos_list = []  # id,start_new,end_new,pos_old
    for i, sample_id in enumerate(id_list):  # id_list
        #     print(sample_id)
        start_end_set = set()
        se_df = startstop_df[startstop_df.index == sample_id].sort_values('length',
                                                                          ascending=False)  # เรียงความยาวมากไปน้อย
        se = pd.unique(se_df['start_end'])  # หาตำแหน่งเริ่มและจบของสายแต่ละสาย
        #         print(se)
        for j in se:
            (a, b) = tuple(map(int, j.split('_')))
            #         print(a,b)
            if (start_end_set) & set(range(a, b + 1)) == set():  # สายไม่ซ้ำกัน
                align_array[i][a - 1:b] = 1
                start_end_set = start_end_set.union(set(range(a, b + 1)))
                choose_pos_list += [[sample_id, a, b, j]]
            else:  # สายซ้ำกัน
                #             print(sample_id)
                intersec_set = (start_end_set) & set(range(a, b + 1))
                new_pos = set(range(a, b + 1)) - intersec_set
                if new_pos != set():
                    #                 print(sample_id)
                    #                 print(sample_id,min(new_pos),max(new_pos),j)
                    choose_pos_list += [[sample_id, min(new_pos), max(new_pos), j]]
                    for k in new_pos:
                        align_array[i][k - 1] = 1

                snps_max_df = drop_data(snps_max_df, sample_id, intersec_set,
                                        j)  # ลบ snps ในเส้นที่สั้นกว่่า ในตำแหน่งที่ชนกัน
                start_end_set = start_end_set.union(set(range(a, b + 1)))

        sum_length += [[sample_id, len(start_end_set)]]
    return sum_length, snps_max_df, align_array, choose_pos_list


# -------------- predict clade ---------------------------------------------------------

# def predict_clade(id_list,snps_max_df,align_array):

#     clade_dict = {}
#     clade_dict['S'] = set(['C8782T', 'T28144C'])
#     clade_dict['L'] = set(['C241', 'C3037', 'C8782', 'G11083', 'C22227','A23063', 'A23403',
#                            'G25563', 'G26144', 'T28144', 'G28882'])  #reference sequence
#     clade_dict['V'] = set(['G11083T', 'G26144T'])
#     clade_dict['G'] = set(['C241T', 'C3037T', 'A23403G'])
#     clade_dict['GH'] = set(['C241T', 'C3037T', 'A23403G', 'G25563T'])
#     clade_dict['GR'] = set(['C241T', 'C3037T', 'A23403G', 'G28882A'])
#     clade_dict['GV'] = set(['C241T', 'C3037T', 'A23403G', 'C22227T'])
#     clade_dict['GRY'] = set(['C241T', 'C3037T', 'A23403G', 'G28882A', 'A23063T',
#                             'T21765-', 'A21766-', 'C21767-', 'A21768-', 'T21769-', 
#                              'G21770-', 'T21991-', 'T21992-', 'A21993-'])
#     clade_list = []
#     clade_result = []

#     # for p in set(snps_maxlength_df['id']):
#     for p in id_list:
#         clade_list += [p]

#         for c in clade_dict.keys():
#             intersec_snp = clade_dict[c].intersection(list(snps_max_df[snps_max_df['id']==p]['change']))
#             len_intersec = len(intersec_snp)
#             clade_list += [len_intersec]
#             clade_list += [intersec_snp]

#         clade_result += [clade_list]
#         clade_list = []

#     clade_df = pd.DataFrame(clade_result, columns=['id','S','snp_S','L','snp_L','V','snp_V','G','snp_G',
#                                                    'GH','snp_GH','GR','snp_GR','GV','snp_GV','GRY','snp_GRY'])
#     #หาความน่่าจะเป็น clade นี้
#     clade_df['S'] = clade_df['S']/2
#     clade_df['V'] = clade_df['V']/2
#     clade_df['G'] = clade_df['G']/3
#     clade_df['GH'] = clade_df['GH']/4
#     clade_df['GR'] = clade_df['GR']/4
#     clade_df['GV'] = clade_df['GV']/4
#     clade_df['GRY'] = clade_df['GRY']/14

#     clade_df['max'] = clade_df[['S','L','V','G','GH','GR','GV','GRY']].max(axis=1)

#     clade_df['result'] = 'O'   #เริ่มต้น ให้ทุก sample ไม่รู้ clade

#     for _,row in clade_df.iterrows():

#         #เช็คตำแหน่ง ใส่ว่าเป็น clade L
#         if all(align_array[_][elem]!=''  for elem in [241,3037,8782,11083,22227,23063,23403,25563,26144,28882,28144]):
#             if row['max']==0:
#                 clade_df.loc[_,'result'] = 'L'

#         #เช็คว่ามี snps ครบหรือไม่ตามตาราง
#         if row['max']==1:
#             if row['GRY'] == row['max']:
#                 clade_df.loc[_,'result'] = 'GRY'
#             elif row['GV'] == row['max']:
#                 clade_df.loc[_,'result'] = 'GV'
#             elif row['GR'] == row['max']:
#                 clade_df.loc[_,'result'] = 'GR'
#             elif row['GH'] == row['max']:
#                 clade_df.loc[_,'result'] = 'GH'
#             elif row['G'] == row['max']:
#                 clade_df.loc[_,'result'] = 'G'
#             elif row['S'] == row['max']:
#                 clade_df.loc[_,'result'] = 'S'
#             elif row['V'] == row['max']:
#                 clade_df.loc[_,'result'] = 'V'
#     return clade_df


def predict_clade(id_list, snps_max_df, align_array, preprocess_dir):
    clade_data = pd.read_csv('clade.tsv', sep='\t')

    clade_dict = clade_data.set_index('Clade').T.to_dict('list')

    # delete nan
    for k in clade_dict.keys():
        clade_dict[k] = [x for x in clade_dict[k] if str(x) != 'nan']

    clade_list = []
    clade_result = []

    # for p in set(snps_maxlength_df['id']):
    for p in id_list:
        clade_list += [p]

        for c in clade_dict.keys():
            intersec_snp = set(clade_dict[c]).intersection(list(snps_max_df[snps_max_df['id'] == p]['change']))
            len_intersec = len(intersec_snp)
            clade_list += [len_intersec]
            clade_list += [intersec_snp]

        clade_result += [clade_list]
        clade_list = []

    # create dataframe
    colums_name = ['id']
    for k in clade_dict.keys():
        colums_name.append(k)
        colums_name.append('snps_{}'.format(k))

    clade_df = pd.DataFrame(clade_result, columns=colums_name)

    # หาความน่่าจะเป็น clade นี้
    for k in clade_dict.keys():
        clade_df[k] = clade_df[k] / len(clade_dict[k])

    clade_df['max'] = clade_df[[k for k in clade_dict.keys()]].max(axis=1)

    clade_df['result'] = 'O'  # เริ่มต้น ให้ทุก sample ไม่รู้ clade

    clade_L_pos = set()
    #     clade_L_pos = list(clade_L_pos.union([int(p[1:])-1 for p in clade_dict['L']]))
    clade_L_pos = list([int(p[1:]) - 1 for p in clade_dict['L']])
    clade_L_neu = list([str(p[0]) for p in clade_dict['L']])

    for _, row in clade_df.iterrows():

        # เช็คตำแหน่ง ใส่ว่าเป็น clade L
        #         if all(align_array[_][elem]!=''  for elem in clade_L_pos):
        if all(align_array[_][p] == n for p, n in zip(clade_L_pos, clade_L_neu)):
            if row['max'] == 0:
                clade_df.loc[_, 'L'] = 1
                clade_df.loc[_, 'max'] = 1
                clade_df.loc[_, 'result'] = 'L'

        # เช็คว่ามี snps ครบหรือไม่ตามตาราง
        if row['max'] == 1:
            if row['GRY'] == row['max']:
                clade_df.loc[_, 'result'] = 'GRY'
            elif row['GH'] == row['max']:
                clade_df.loc[_, 'result'] = 'GH'
            elif row['GR'] == row['max']:
                clade_df.loc[_, 'result'] = 'GR'
            elif row['GV'] == row['max']:
                clade_df.loc[_, 'result'] = 'GV'
            elif row['G'] == row['max']:
                clade_df.loc[_, 'result'] = 'G'
            elif row['S'] == row['max']:
                clade_df.loc[_, 'result'] = 'S'
            elif row['V'] == row['max']:
                clade_df.loc[_, 'result'] = 'V'
    return clade_df


# -------------- check protein ---------------------------------------------------------

from Bio.Seq import Seq


# check ว่า protein เปลี่ยนรึเปล่่า
def check_protein(id_number, p, gene_name, df, position_gene_df, seq_ref_list):
    gene_start = position_gene_df[position_gene_df['gene'] == gene_name]['start'].values
    gene_start = gene_start - 1
    #     mod = (p-gene_start-1)%3
    #     sample_df = df.loc[df.index==id_number,:].values[0]

    #     if mod == 1:
    #         if (sample_df[p-1]!='-') & (sample_df[p]!='-') & (sample_df[p+1]!='-'):
    #             s = sample_df[p-1]+sample_df[p]+sample_df[p+1]
    #             r = seq_ref_list[p-1]+seq_ref_list[p]+seq_ref_list[p+1]
    #             check = (Seq(s).translate()==Seq(r).translate())
    #             s_p = str(Seq(s).translate())
    #             r_p = str(Seq(r).translate())
    #         else:
    #             check = 'unknown'
    #             s_p = '-'
    #             r_p = '-'
    #     elif mod == 2:
    #         if (sample_df[p-2]!='-') & (sample_df[p-1]!='-') & (sample_df[p]!='-'):
    #             s = sample_df[p-2]+sample_df[p-1]+sample_df[p]
    #             r = seq_ref_list[p-2]+seq_ref_list[p-1]+seq_ref_list[p]
    #             check = (Seq(s).translate()==Seq(r).translate())
    #             s_p = str(Seq(s).translate())
    #             r_p = str(Seq(r).translate())
    #         else:
    #             check = 'unknown'
    #             s_p = '-'
    #             r_p = '-'
    #     elif mod == 0:
    #         if (sample_df[p]!='-') & (sample_df[p+1]!='-') & (sample_df[p+2]!='-'):
    #             s = sample_df[p]+sample_df[p+1]+sample_df[p+2]
    #             r = seq_ref_list[p]+seq_ref_list[p+1]+seq_ref_list[p+2]
    #             check = (Seq(s).translate()==Seq(r).translate())
    #             s_p = str(Seq(s).translate())
    #             r_p = str(Seq(r).translate())
    #         else:
    #             check = 'unknown'
    #             s_p = '-'
    #             r_p = '-'
    #     else: print('error')

    #     return check, r_p, s_p, int((p-gene_start-1)/3)+1

    mod = (p - gene_start - 1) % 3
    sample_seq = df.loc[df.index == id_number, :].values[0]
    codon = None
    if mod == 1:
        x, y, z = p - 1, p, p + 1
    elif mod == 2:
        x, y, z = p - 2, p - 1, p
    elif mod == 0:
        x, y, z = p, p + 1, p + 2

    s = sample_seq[x] + sample_seq[y] + sample_seq[z]
    r = seq_ref_list[x] + seq_ref_list[y] + seq_ref_list[z]
    r_p = str(Seq(r).translate())

    if (sample_seq[x] == '-' and sample_seq[y] == '-' and sample_seq[z] == '-') or (
            sample_seq[x] != '-' and sample_seq[y] != '-' and sample_seq[z] != '-'):
        check = (Seq(s).translate(gap='-') == Seq(r).translate())
        s_p = str(Seq(s).translate(gap='-'))
        # if '---' then '-', if there exist alphabet not in {ATCG}, then 'X'
    else:
        check = False
        s_p = 'X'

    # remove duplicate results will be removed in the later step
    return check, r_p, s_p, int((p - gene_start - 1) / 3) + 1


def change_protein_table(change_df, position_gene_df, align_array_df, seq_ref_list):
    # ตารางการเปลี่ยนกรดอะมิโน ในโปรตีนต่างๆ
    change_protein = []
    gene_dict = {}

    for i, p in enumerate(change_df['position']):

        gene = 'None'
        id_number = change_df.loc[i, 'id']
        #     print(id_number)
        for _, row in position_gene_df.iterrows():
            if p >= row['start'] and p <= row['end']:
                gene = row['gene']
            #     print(i,p,gene)
        if gene != 'None':
            #         print(i,p, gene)
            check, ref_protein, sample_protein, pos_protein = check_protein(id_number, p, gene, align_array_df,
                                                                            position_gene_df, seq_ref_list)
            #         print(check, ref_protein, sample_protein, pos_protein)
            change_protein += [[id_number, p, ref_protein, pos_protein, sample_protein,
                                ref_protein + str(pos_protein) + sample_protein, gene, check]]
        gene_dict[p] = gene
        #             print(p,check)

    change_protein_df = pd.DataFrame(change_protein,
                                     columns=['id', 'pos_neucleotide', 'query_protein', 'position_protein',
                                              'sbjct_protein', 'change_protein', 'gene', 'check'])
    change_protein_df['pos_neucleotide'] += 1

    return change_protein_df


# ------------------------text in dendro ----------------------------------
def get_mutation_text(m_list):
    m_list = [m for m in m_list if (m[-1] != 'N') and (m[-1] != '-') and ('*' not in m)]
    #     print(m_list[0][-1])

    if len(m_list) == 0:
        return '-'
    elif len(m_list) <= 3:
        return ', '.join(m_list)
    else:
        m_text = ''
        for i, m in enumerate(m_list):
            m_text += m + ', '
            if (i % 3) == 2:
                m_text += '</br> '
        if len(m_list) % 3 == 0:
            return m_text[:-8]
        else:
            return m_text[:-2]


# ------------------------------------------------------------------
def find_align_array_pro(align_array_df, position_gene_df, change_df):
    id_list = list(align_array_df.index)

    # create dict_protein
    position_gene_df['start-end'] = [(a, b) for a, b in zip(position_gene_df['start'], position_gene_df['end'])]
    position_gene_df = position_gene_df[['gene', 'start-end']]
    dict_protein = dict(zip(position_gene_df['gene'], position_gene_df['start-end']))

    align_array_pro = align_array_df.copy()
    A = align_array_pro.to_numpy()

    for name, r in list(dict_protein.items()):
        i, j = r
        for k in range(len(id_list)):

            neu = A[k][i:j]
            #         print(name,i,j,len(neu))
            pos = np.where(neu == '-')
            list_pos = list(pos[0])  # ตำแหน่งที่่มี deletion ทั้งหมด
            #         print(list_pos)

            del_list = []
            start_del_list = []
            num_del_list = []
            for p in list_pos:
                if (p not in del_list) and (np.isin([p, p + 1, p + 2], list_pos).sum() == 3):
                    if del_list == []:
                        start_del_list += [p]  # ตัวแรก
                    if del_list != [] and (p - del_list[-1]) != 1:
                        start_del_list += [p]  # ตัวที่ไม่ต่อกัน
                        num_del_list += [len(del_list)]  # ความถี่สะสม

                    del_list += [p, p + 1, p + 2]
            num_del_list += [len(del_list)]
            #             if start_del_list!=[]:
            #                 print(k,name,start_del_list,num_del_list)

            # เปลี่ยนใน array
            for m in range(len(start_del_list)):
                pos_start = start_del_list[m] + dict_protein[name][0]
                pos_end = start_del_list[m] + num_del_list[m] + dict_protein[name][0]

                if start_del_list[m] % 3 == 2:
                    #                     print('change0:{}={}'.format(name,pos_start))
                    align_array_pro.iloc[k][pos_start] = align_array_pro.iloc[k][pos_end]
                    align_array_pro.iloc[k][pos_end] = '-'

                    temp_list = [[id_list[k], pos_end, 'None']]
                    #                     print(temp_list)
                    temp_df = pd.DataFrame(temp_list, columns=(['id', 'position', 'check']))
                    change_df = pd.concat([change_df, temp_df])

                elif start_del_list[m] % 3 == 1:
                    #                     print('change1:{}={}'.format(name,pos_start))
                    align_array_pro.iloc[k][pos_start] = align_array_pro.iloc[k][pos_end]
                    align_array_pro.iloc[k][pos_start + 1] = align_array_pro.iloc[k][pos_end + 1]
                    align_array_pro.iloc[k][pos_end] = '-'
                    align_array_pro.iloc[k][pos_end + 1] = '-'

                    temp_list = [[id_list[k], pos_end, 'None'], [id_list[k], pos_end + 1, 'None']]
                    #                     print(temp_list)
                    temp_df = pd.DataFrame(temp_list, columns=(['id', 'position', 'check']))
                    change_df = pd.concat([change_df, temp_df])

    return align_array_pro, change_df
