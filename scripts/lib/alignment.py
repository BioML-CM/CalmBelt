from Bio.Blast import NCBIXML
import pandas as pd
import pickle



# def hsp_to_snps(hsp,person):
#     q = hsp.query  #reference
#     ms = hsp.match                
#     s = hsp.sbjct  #patient
#     start = hsp.query_start             
#     end = hsp.query_end
    
#     snps_list = []
    
#     for i,m in enumerate(ms):
#         if m!='|':
#             pos = i+start      
# #             snps_list += ['{}{}{}'.format(q[i],pos,s[i])]
#             snps_list +=[[person, q[i], pos, s[i], "{}_{}".format(start,end), end-start+1]]
#     return snps_list



# def start_to_end(hsp,person):
#     start = hsp.query_start             
#     end = hsp.query_end
#     return [person, "{}_{}".format(start,end), end-start+1]


# def start_end_list_func(url):
#     result = []
#     start_end = []
#     snps_list = []
#     start_end_list = []
#     match = []


#     for record in NCBIXML.parse(open(url)): 
#         print(record)
#         if record.alignments:         
#             for i, align in enumerate(record.alignments): 
#                 person = str(align).split('|')[1].split('_')[2]     
#                 match += [person]


#                 for hsp in align.hsps:                
#                     snps_list += hsp_to_snps(hsp,person)
#                     start_end_list += [start_to_end(hsp,person)]


#                 result += snps_list
#                 start_end += start_end_list
#                 snps_list = []           
#                 start_end_list = []
                
#     return result,start_end



# def drop_data(snps_max_df,sample_id,intersec_set,start_end):
#     snps_id = snps_max_df[snps_max_df['id']==sample_id]    
#     position = list(snps_id[snps_id['start_end']==start_end]['position'])
#     for pos in position:
#         if pos in intersec_set:
#             indexNames = snps_max_df[(snps_max_df['start_end'] == start_end) & (snps_max_df['position'] == pos)].index
#             snps_max_df.drop(indexNames , inplace=True)
#     return snps_max_df



# def drop_query(snps_max_df,sample_id,start_end):
#     snps_id = snps_max_df[snps_max_df['id']==sample_id]    
# #     position = list(snps_id[snps_id['start_end']==start_end]['position'])
# #     for pos in position:
# #         if pos in intersec_set:
#     indexNames = snps_max_df[(snps_max_df['start_end'] == start_end) & (snps_max_df['query'] == '-')].index
#     snps_max_df.drop(indexNames , inplace=True)
#     return snps_max_df


# def deletion(start_end_del,snps_df,sample_id):
#     for sed in start_end_del:
#         new_pos = []
#         count = 1
#         i=0

#         pos_del = sorted(list(snps_df[(snps_df['id']==sample_id)&(snps_df['query']=='-')&(snps_df['start_end']==sed)]['position']))
#         snp_del = snps_df[(snps_df['id']==sample_id)&(snps_df['start_end']==sed)].sort_values(['position']).copy()       
#         pos_all = sorted(list(snp_del['position']))



#         for a in pos_all:   
#             old_pos = snp_del[snp_del['position']==a]['position'].values[0]
#             if a<=pos_del[i]: #ไม่เปลี่ยน
#                 new_pos += [old_pos]
#             else:
#                 if pos_del[i]==pos_del[-1]:     #a>pos_del[i] และ pos_del[i]เป็นตัวสุดท้าย
#                     new_pos += [old_pos-count]
#                 else:
#                     if a<pos_del[i+1]:              #a>pos_del[i] และ pos_del[i]ไมใช่ตัวสุดท้าย
#                         new_pos += [old_pos-count]
#                     else : 
#                         count += 1
#                         i+= 1       
#                         new_pos += [old_pos-count]

#         snps_df.loc[(snps_df['id']==sample_id)&(snps_df['start_end']==sed),'position'] = new_pos        
#     #         print(sed,new_pos)
#         snps_df = drop_query(snps_df,sample_id,sed)
#     return snps_df


# def length_max_snps(id_list,startstop_df,snps_max_df,align_array):
#     sum_length = []
#     for i,sample_id in enumerate(id_list): #id_list
#         start_end_set = set()
#         se_df = startstop_df[startstop_df['id'] ==sample_id].sort_values('length', ascending=False)
#         se = pd.unique(se_df['start_end'])   #list(set(se_df['start_end']))
#         for j in se: 
#             (a, b) = tuple(map(int, j.split('_')))
#             if (start_end_set) & set(range(a,b+1)) == set(): #สายไม่ซ้ำกัน
#                 align_array[i][a-1:b] =1
#                 start_end_set = start_end_set.union(set(range(a,b+1)))
#             else: #สายซ้ำกัน
#                 intersec_set = (start_end_set) & set(range(a,b+1))
#                 new_pos = set(range(a,b+1))-intersec_set
#                 print(new_pos)
#                 for k in new_pos:
#                     align_array[i][k-1] =1
#                 print('id: {} intersec!!! position {}'.format(sample_id,j))
#                 snps_max_df = drop_data(snps_max_df,sample_id,intersec_set,j)
#                 start_end_set = start_end_set.union(set(range(a,b+1)))
#                 print(len(start_end_set))
#         sum_length += [[sample_id, len(start_end_set)]]
    
#     return sum_length,snps_max_df



# from Bio.Seq import Seq
# #check ว่า protein เปลี่ยนรึเปล่่า
# def check_protein(id_number, p, gene_name, df, position_gene_df, seq_ref_list):
#     gene_start = position_gene_df[position_gene_df['gene']==gene_name]['start'].values
#     gene_start = gene_start-1
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

# TODO remove N from m_list (for DNA-level)
# - DNA, show only SNPS
# - Protein, show amino changes and '-'
def get_mutation_text(m_list):
    m_list = [m for m in m_list if (m[-1]!='N') and (m[-1]!='-') and ('*' not in m)]
#     print(m_list[0][-1])
    
    if len(m_list) == 0 :  
        return '-'
    elif len(m_list)<=3:
        return ', '.join(m_list) 
    else:
        m_text = ''
        for i,m in enumerate(m_list):            
            m_text += m+', '
            if (i%3)==2:
                m_text += '</br> '
        if len(m_list)%3==0:
            return m_text[:-8]
        else: 
            return m_text[:-2]