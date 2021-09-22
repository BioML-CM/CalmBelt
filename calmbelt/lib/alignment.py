from Bio.Blast import NCBIXML
import pandas as pd
import pickle


def get_mutation_text(m_list):
    m_list = [m for m in m_list if (m[-1]!='N') and (m[-1]!='-') and ('*' not in m)]
    
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