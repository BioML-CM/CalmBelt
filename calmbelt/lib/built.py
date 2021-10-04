import pandas as pd
import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt

import plotly.express as px
import plotly.graph_objects as go

from lib import alignment

from plotly.figure_factory import create_dendrogram
from scipy.spatial.distance import pdist

import plotly.figure_factory as ff
import itertools

from dateutil.relativedelta import relativedelta


# ----------------------------------------------------------------
def get_tick_by_time(scale, min_date, max_date):
    current = min_date
    result = [current + relativedelta(months=-1)]

    while current <= max_date:
        result.append(current)
        current += relativedelta(months=1)
    ticktext = [i.strftime('%b-%Y') for i in result]

    tickvals = [scale * i for i in range(len(ticktext))]

    return tickvals, ticktext


# ----------------------------------------------------------------

def find_insertion_position(id_insertion, snps_insertion_df):
    insertion_list = []
    for sample in id_insertion:
        se_list = list(set(snps_insertion_df[snps_insertion_df['id'] == sample]['start_end']))

        for se in se_list:  # แยกคิดแต่ละสาย
            list_temp = sorted(list(
                snps_insertion_df[(snps_insertion_df['id'] == sample) & (snps_insertion_df['start_end'] == se)][
                    'position']))
            list_temp += [0]

            count_insertion = []  # จำนวนตำแหน่งที่่อยู่ติดกัน
            pos_temp = []
            sum_count = 1
            for i in range(len(list_temp) - 1):
                if list_temp[i] == list_temp[i + 1] - 1:
                    sum_count += 1
                else:
                    count_insertion += [sum_count]  # จำนวน
                    pos_temp += [list_temp[i]]  # ตำแหน่ง
                    sum_count = 1

            # ตำแหน่งตาม ref
            pos_insertion = [i - j for i, j in zip(pos_temp, np.cumsum(count_insertion))]
            insertion_list += [[sample, i, j, se] for i, j in zip(pos_insertion, count_insertion)]
    return insertion_list


# ----------------------------------------------------------------------


def fig_dendro_built_iqtree(X_dendro, X_id_list, clade_drop_dendro, snps_max_df, change_protein_df):
    n = len(X_id_list)

    fig_dendro = create_dendrogram(X_dendro, orientation='right', labels=X_id_list,
                                   distfun=lambda x: pdist(x, metric='hamming'))

    labels_list = fig_dendro.layout['yaxis']['ticktext']  # label เรียงตามลำดับจากล่างขึ้นบน
    # สร้าง heatmap
    m = np.zeros([n, n])
    vals = pdist(X_dendro, metric='hamming')

    xs, ys = np.triu_indices(n, k=1)
    m[xs, ys] = vals
    m[ys, xs] = vals
    m[np.diag_indices(n)] = 0
    m = m * 29903  # จำนวน dna
    m = np.round(m, 0).astype(int)

    new_index = [X_id_list.index(k) for k in list(labels_list)]

    m = m[new_index]  # row=yaxis
    m = m[:, new_index[::-1]]  # column=xaxis

    fig_heatmap = ff.create_annotated_heatmap(m, x=list(labels_list[::-1]), y=list(labels_list), colorscale='cividis',
                                              hoverinfo='none')

    heatmap_df = pd.DataFrame(m[::-1], index=list(labels_list[::-1]), columns=list(labels_list[::-1]))

    # Make text size smaller
    for i in range(len(fig_heatmap.layout.annotations)):
        fig_heatmap.layout.annotations[i].font.size = 9

    h = max((int(n / 20) + 1) * 300, 600)
    fig_heatmap['layout'].update({'template': 'plotly_white', 'height': h, 'width': h,
                                  'xaxis': {'title': '', 'showticklabels': True, 'ticks': ''},
                                  'yaxis': {'title': '', 'showticklabels': True, 'ticks': ''},
                                  })

    # -----------------------------

    position_dendro = []

    size = len(fig_dendro.data)

    min_x = 0
    for i in range(size):
        x = fig_dendro.data[i]['x']
        y = fig_dendro.data[i]['y']
        fig_dendro.data[i]['marker']['color'] = 'grey'
        fig_dendro.data[i]['marker']['size'] = 0.2
        fig_dendro.data[i]['showlegend'] = False

        min_x = min(min(x), min_x)

    position_dendro_df = pd.DataFrame(labels_list, columns=['id'])
    position_dendro_df = position_dendro_df.sort_values(['id'])
    position_dendro_df = position_dendro_df.reset_index()

    position_dendro_df = pd.concat(
        [position_dendro_df, clade_drop_dendro[['year', 'month', 'day', 'group', 'lineage', 'result']]], axis=1)

    # หาเวลาที่น้อยที่่สุดของข้อมูล
    position_dendro_df['date'] = pd.to_datetime(position_dendro_df[['year', 'month', 'day']])

    df = position_dendro_df[position_dendro_df['id'] != 'Reference-30-12-2019']
    min_year = pd.to_datetime(df[['year', 'month', 'day']]).min().strftime('%Y')
    min_month = pd.to_datetime(df[['year', 'month', 'day']]).min().strftime('%m')
    min_day = pd.to_datetime(df[['year', 'month', 'day']]).min().strftime('%d')

    df['date'] = pd.to_datetime(df[['year', 'month', 'day']])
    min_date = df['date'].min()
    max_date = df['date'].max()

    position_dendro_df.loc[position_dendro_df['id'] == 'Reference-30-12-2019', 'year'] = min_year
    position_dendro_df.loc[position_dendro_df['id'] == 'Reference-30-12-2019', 'month'] = min_month
    position_dendro_df.loc[position_dendro_df['id'] == 'Reference-30-12-2019', 'day'] = min_day
    position_dendro_df.loc[position_dendro_df['id'] == 'Reference-30-12-2019', 'date'] = min_date

    position_dendro_df['diff'] = position_dendro_df['date'] - min_date

    # วาดเส้นแนวนอน
    x_new = []

    ratio = min(max(3, (int(n / 5) + 2)), 5)
    scale = np.abs(min_x) / ratio
    #     scale = 0.0005
    for _, row in position_dendro_df.iterrows():
        x_new += [((int(row['diff'].days) / 30) * scale) + scale]


    position_dendro_df['xnew'] = x_new
    position_dendro_df = position_dendro_df.sort_values(['index'])


    indications = sorted(list(set(clade_drop_dendro['group'])))

    clade_color = {0: 'lightcoral',
                   1: 'darkcyan',
                   2: 'goldenrod',
                   3: 'olive',
                   4: 'turquoise',
                   5: 'deepskyblue',
                   6: 'slategrey',
                   7: 'yellowgreen',
                   8: 'hotpink',
                   9: 'indianred',
                   10: 'sandybrown',
                   11: 'darkorange',
                   12: 'plum',
                   13: 'salmon',
                   14: 'steelblue',
                   15: 'cornflowerblue',
                   16: 'slateblue',
                   17: 'lightpink',
                   18: 'lightgreen',
                   19: 'dodgerblue',
                   20: 'darkorchid', }

    group_dict = dict(zip(indications, range(len(indications))))
    group_color_dict = {}
    for k, v in group_dict.items():
        group_color_dict[k] = clade_color[v]

    color_list = []
    label_list = []
    for i in X_id_list:
        c = clade_drop_dendro[clade_drop_dendro['id'] == i]['group'].values[0]
        label_list += ['{}'.format(i)]
        color_list += [clade_color[group_dict[c]]]

    color_label = zip(label_list, color_list)
    dict_color_label = dict(color_label)

    color_list_sample = [dict_color_label[str(i)] for i in labels_list]
    position_dendro_df['color'] = color_list_sample

    y_new = []

    for i in range(0, len(labels_list)):
        color_line = dict_color_label[str(labels_list[i])]
        xnew = position_dendro_df[position_dendro_df['index'] == i]['xnew'].values[0]
        fig_dendro.add_trace(
            go.Scatter(x=[0, xnew], y=[5 + i * 10, 5 + i * 10], mode='lines', hoverinfo='none', opacity=0.5,
                       showlegend=False,
                       marker=dict(color=color_line, line=dict(width=0.01))
                       ))
        y_new += [5 + i * 10]

    # วาดจุด
    position_dendro_df['ynew'] = y_new

    change_protein_df = change_protein_df[change_protein_df['id'].isin(X_id_list)][
        ['id', 'change_protein', 'gene', 'check']]

    for i in range(n):
        x = position_dendro_df[position_dendro_df['id'] == labels_list[i]]['xnew']  # [position_dendro_df['xnew'][i]]
        y = position_dendro_df[position_dendro_df['id'] == labels_list[i]]['ynew']  # [position_dendro_df['ynew'][i]]
        id_sample = str(labels_list[i])

        snps_max_df = snps_max_df.sort_values('position')
        change_neucleotide_list = list(snps_max_df[snps_max_df['id'] == labels_list[i]]['change'])
        change_neucleotide = alignment.get_mutation_text(change_neucleotide_list)

        change_protein_label = change_protein_df[
            (change_protein_df['id'] == labels_list[i]) & (change_protein_df['check'] == False)]
        if len(list(change_protein_label['change_protein'])) == 0:
            change_protein = '-'
        else:
            change_protein = ''
            gene = list(set(change_protein_label['gene']))
            for g in ['S']:  # gene:
                list_gene = list(set(change_protein_label[change_protein_label['gene'] == g]['change_protein']))

                change_protein += f'<i>{g}</i>' + ' : ' + alignment.get_mutation_text(list_gene) + '</br> '

        color_marker = dict_color_label[str(id_sample)]
        clade = position_dendro_df[position_dendro_df['id'] == labels_list[i]]['result'].values[0]
        lin = position_dendro_df[position_dendro_df['id'] == labels_list[i]]['lineage'].values[0]
        text_lin = id_sample + ' (' + lin + ')'
        g = position_dendro_df[position_dendro_df['id'] == labels_list[i]]['group'].values[0]

        fig_dendro.add_trace(go.Scatter(x=x, y=y, text=text_lin,
                                        mode='markers+text', name='', textposition="middle right",
                                        textfont=dict(size=8, ),
                                        hovertemplate='<b>Id</b> : {}'.format(id_sample) +
                                                      '<br><b>Group</b> : {}'.format(g) +
                                                      '<br><b>Lineage</b> : {}'.format(lin) +
                                                      '<br><b>Clade</b> : {}'.format(clade) +
                                                      '<br><b>Mutation</b>:<br> {}'.format(change_neucleotide) +
                                                      '<br><br><b>Amino acid change (Spike) </b>:<br> {}'.format(
                                                          change_protein),
                                        hoverlabel=dict(bgcolor='#F2F2F3'),
                                        marker=dict(color=color_marker, size=8, line_width=1),
                                        legendgroup=g, showlegend=False))

    for k, v in group_dict.items():
        fig_dendro.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                                        marker=dict(size=8, color=clade_color[v], line_width=1),
                                        legendgroup=k, showlegend=True, name=k))

    tickvals, ticktext = get_tick_by_time(scale, min_date, max_date)
    fig_dendro.update_layout(
        xaxis=dict(
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext,
            tickangle=-45,
        ),
        margin=dict(l=0, r=0, b=0, t=50, pad=4)
    )
    fig_dendro.update_xaxes(showgrid=True)

    h = max(500, (int(n / 30) + 1) * 300)

    fig_dendro['layout'].update({'template': 'plotly_white',
                                 'height': h,
                                 'xaxis': {'title': '', 'showticklabels': True, 'ticks': ''},
                                 'yaxis': {'title': '', 'showticklabels': False, 'ticks': ''},
                                 })
    fig_dendro.update_layout(legend_font_size=10)
    fig_dendro.update_xaxes(tickfont_size=10)

    fig_dendro.update_xaxes(range=(1.5 * min_x, position_dendro_df['xnew'].max() + np.abs(min_x) * 1.5))
    fig_dendro.update_yaxes(range=(-20, position_dendro_df['ynew'].max() + 25))
    fig_dendro.update_layout(showlegend=True)
    fig_dendro['layout'].update({'legend': {'traceorder': 'normal'}})

    return fig_dendro, fig_heatmap, position_dendro_df, group_color_dict, heatmap_df


# -------------------------------------------------------
import toytree  # a tree plotting library
import toyplot  # a general plotting library
import toyplot.png
from Bio import SeqIO
import glob
import os
import subprocess


def built_iqtree(group_dict, path):

    record_list = []
    for name in glob.glob(os.path.join(path, '*_*_*.fasta')):
        for record in SeqIO.parse(name, "fasta"):
            record.id = '_'.join(name.split('/')[-1].replace('.fasta', '').split('_')[2:5])
            record_list += [record]

    with open(os.path.join(path, 'concat.fasta'), "w") as output_handle:
        SeqIO.write(record_list, output_handle, "fasta")

    command_mafft_list = ['mafft', '--retree', '2', f'{os.path.join(path, "concat.fasta")}', '>',
                          f'{os.path.join(path, "concat.mafft")}']
    mafft_output = subprocess.check_output(' '.join(command_mafft_list), shell=True)

    command_iqtree_list = ['iqtree2', '-s', f'{os.path.join(path, "concat.mafft")}', '-m', 'GTR+G']
    iqtree_output = subprocess.check_output(' '.join(command_iqtree_list), shell=True)

    newick_fname = os.path.join(path, 'concat.mafft.treefile')
    rtre = toytree.tree(newick_fname, tree_format=1)

    colorlist = [group_dict[i.split('_')[-1]] for i in rtre.get_tip_labels()]
    canvas, axes, mark = rtre.draw(tip_labels_colors=colorlist, scalebar=True)

    # save png
    toyplot.png.render(canvas, f'{os.path.join(path, "iqtree.png")}', scale=3)


# -------------------------------------------------------

def built_timetree(group_dict, path):

    aln_fname = os.path.join(path, "concat.mafft")
    tree_fname = os.path.join(path, "concat.mafft.treefile")
    new_newick_fname = os.path.join(path, "concat.mafft.timetree.nwk")

    sample_date_dict = {}
    with open(aln_fname) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == 'Reference':
                sy = int(date(2019, 12, 30).strftime("%Y"))
                sd = (float(date(2019, 12, 30).strftime("%j")) - 1) / 365
                sdate = sy + sd
            else:
                sid, sdate, group = record.id.split('_')
                d, m, y = sdate.split('-')
                sy = int(date(int(y), int(m), int(d)).strftime("%Y"))
                sd = (float(date(int(y), int(m), int(d)).strftime("%j")) - 1) / 365
                sdate = sy + sd

            sample_date_dict[record.id] = sdate

    tt = TreeTime(tree=tree_fname, aln=aln_fname, dates=sample_date_dict)
    tt.run(root='Reference')

    Phylo.write(tt.tree, new_newick_fname, "newick")

    rtre = toytree.tree(new_newick_fname, tree_format=1)

    colorlist = [group_dict[i.split('_')[-1]] for i in rtre.get_tip_labels()]
    canvas, axes, mark = rtre.draw(tip_labels_colors=colorlist, scalebar=True)

    # save png
    toyplot.png.render(canvas, f'{os.path.join(path, "timetree.png")}', scale=3)

