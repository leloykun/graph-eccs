import streamlit as st
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np

st.title("Tree Eccentricity Visualizer")

version = st.sidebar.selectbox(
    "Version",
    ("Standard Unlabaled Tree", "Eccentricity Sequence Isomorphism"),
    index=1
)
MIN_N = 3
MAX_N = 21 if version == "Standard Unlabaled Tree" else 30

@st.cache(suppress_st_warning=True, allow_output_mutation=True)
def load_prufer_codes(n):
  prufer_codes_ = []
  if version == "Standard Unlabaled Tree":
    filename = "trees/tree_prufer_{}.txt".format(n)
  else:
    filename = "trees/tree_seq_iso_{}.txt".format(n)
  with open(filename, 'r') as f:
    N, M = map(int, f.readline().split())
    for _ in range(M):
      code = list(map(int, f.readline().split()))
      prufer_codes_.append(code)
  return prufer_codes_

@st.cache(suppress_st_warning=True, allow_output_mutation=True)
def load_ecc_data(n):
  bis_ = []
  seq_ = []
  perm_ = []
  if version == "Standard Unlabaled Tree":
    filename = "trees/ecc_{}.txt".format(n)
  else:
    filename = "trees/ecc_seq_{}.txt".format(n)
  with open(filename, 'r') as f:
    N, M = map(int, f.readline().split())
    for i in range(M):
      line = f.readline().split()
      bis_.append(int(line[0]))
      seq_.append(line[1:-1])
      perm_.append(int(line[-1]))
  return bis_, seq_, perm_

@st.cache(suppress_st_warning=True, allow_output_mutation=True)
def load_ecc2data(n):
  prufer_codes_ = load_prufer_codes(n)
  bis_, ecc_seq_, perm_ = load_ecc_data(n)
  print(len(perm_), len(ecc_seq_), len(prufer_codes_))
  print(perm_)

  bi2data_ = {}
  for i in range(len(perm_)):
    data = [prufer_codes_[perm_[i]], ecc_seq_[i]]
    if bis_[i] in bi2data_:
      bi2data_[bis_[i]].append(data)
    else:
      bi2data_[bis_[i]] = [data]
  return bi2data_

tree_size = st.sidebar.slider("Pick Tree Size", MIN_N, MAX_N)

prufer_codes = load_prufer_codes(tree_size)
bis, ecc_seq, perm = load_ecc_data(tree_size)
bi2data = load_ecc2data(tree_size)

min_bi = min(bis)
max_bi = max(bis)

if min_bi == max_bi:
    st.sidebar.write("Bounding Index: {}".format(min_bi))
    bi = min_bi
else:
    bi = st.sidebar.slider("Pick Bounding Index", min_bi, max_bi)

all_accounted_for = np.all(np.isin(np.unique(bis), np.arange(min_bi, max_bi+1)))
if all_accounted_for:
    st.sidebar.write("All bounding indices are accounted for.")
else:
    isin = np.isin(np.arange(min_bi, max_bi+1), np.unique(bis))
    st.sidebar.write("Missing:")
    st.sidebar.write(isin.astype(int))
st.sidebar.write("---")

num_trees = len(bi2data[bi])
if num_trees == 1:
  st.sidebar.write("## There is 1 tree of size {}".format(tree_size))
else:
  st.sidebar.write("## There are {} trees of size {}".format(len(prufer_codes), tree_size))
st.sidebar.write("## Of which, {} has Bounding Index {}".format(num_trees, bi))
st.sidebar.write("---")

st.sidebar.write("Samples")
sample_idx_from = st.sidebar.number_input(
  "From (0 -> {})".format(num_trees-1),
  min_value=0,
  max_value=num_trees-1,
  step=5
)
sample_idx_to = st.sidebar.number_input(
  "To ({} -> {})".format(sample_idx_from, num_trees-1),
  min_value=sample_idx_from,
  max_value=num_trees-1,
  value=min(sample_idx_from+5, num_trees-1),
  step=5
)

st.sidebar.button("Resample")
st.sidebar.write("---")
layout = st.sidebar.selectbox(
    "Layout",
    ("Kamada Kawai", "Spring", "Planar", "Shell"),
    index=3
)
include_label = st.sidebar.checkbox("Include node eccentricity", True)

for i in range(sample_idx_from, sample_idx_to+1):
    code, seq = bi2data[bi][i]
    T = nx.from_prufer_sequence(code)
    fig, ax = plt.subplots()

    non_leaves = []
    T_ = T.copy()
    T_.remove_nodes_from([x for x in T if T.degree(x) == 1])
    for cc in nx.connected_components(T_):
        for s in cc:
            if T_.degree(s) == 1:
                break
        for u in nx.topological_sort(nx.dfs_tree(T_, source=s)):
            non_leaves.append(u)

    leaves = []
    for u in non_leaves:
        for v in T.neighbors(u):
            if v not in non_leaves + leaves:
                leaves.append(v)

    #print("non leaves: ", non_leaves)
    #print("leaves: ", leaves)

    eccs = nx.eccentricity(T)
    ecc_sum = 0
    eccs_list = []
    node_color = []
    for u in eccs:
        eccs_list.append(eccs[u])
        ecc_sum += eccs[u]
    for u in eccs:
        if eccs[u] == min(eccs_list):
            node_color.append('#009FB7')
        elif eccs[u]*len(eccs) < ecc_sum:
            node_color.append('#80A4ED')
        else:
            node_color.append('#F63366')
    eccs_list = sorted(eccs_list, reverse=True)
    bi_ = -1
    for j in range(len(eccs_list)):
        if eccs_list[j]*len(eccs_list) < ecc_sum:
            bi = j
            break
    #print(eccs_list, bi_)

    if layout == "Kamada Kawai":
        pos = nx.kamada_kawai_layout(T)
    elif layout == "Spring":
        pos = nx.spring_layout(T)
    elif layout == "Planar":
        pos = nx.planar_layout(T)
    elif layout == "Shell":
        pos = nx.shell_layout(T, nlist=[non_leaves, leaves])
        pos = nx.spring_layout(T, pos=pos, fixed=non_leaves)

    nx.draw(T, pos=pos, ax=ax, node_color=node_color, with_labels=False)
    if include_label:
        nx.draw_networkx_labels(T, pos=pos, labels=eccs)

    ax.set_title(
        "Prufer Code: [{}]\nAverage Eccentricity: {:.2f}\nEccentricity Sequence: {}".format(
            " ".join(map(str, code)),
            ecc_sum/tree_size,
            " ".join(seq)
        )
    )
    st.pyplot(fig)
plt.clf()
    #dot = nx.nx_pydot.to_pydot(g)
    #st.graphviz_chart(dot.to_string
