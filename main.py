import streamlit as st
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np


@st.cache(suppress_st_warning=True, allow_output_mutation=True)
def load_prufer_codes(n):
  prufer_codes_ = []
  with open("trees/tree_prufer_{}.txt".format(n), 'r') as f:
    N, M = map(int, f.readline().split())
    for _ in range(M):
      code = list(map(int, f.readline().split()))
      prufer_codes_.append(code)
  return prufer_codes_

@st.cache(suppress_st_warning=True, allow_output_mutation=True)
def load_bis(n):
  bis_ = []
  with open("trees/ecc_{}.txt".format(n), 'r') as f:
    N, M = map(int, f.readline().split())
    for i in range(M):
      bi = int(f.readline().split()[-1])
      bis_.append(bi)
  return bis_

@st.cache(suppress_st_warning=True, allow_output_mutation=True)
def load_ecc2graph(n):
  prufer_codes_ = load_prufer_codes(n)
  bis_    = load_bis(n)
  assert(len(prufer_codes_) == len(bis_))

  bi2prufer_codes_ = {}
  for i in range(len(prufer_codes_)):
    if bis_[i] in bi2prufer_codes_:
      bi2prufer_codes_[bis_[i]].append(prufer_codes_[i])
    else:
      bi2prufer_codes_[bis_[i]] = [prufer_codes_[i]]
  return bi2prufer_codes_


st.title("Tree Eccentricity Visualizer")

tree_size = st.sidebar.slider("Pick Tree Size", 3, 21)

prufer_codes = load_prufer_codes(tree_size)
bis = load_bis(tree_size)
bi2prufer_codes = load_ecc2graph(tree_size)

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

num_trees = len(bi2prufer_codes[bi])
if num_trees == 1:
  st.sidebar.write("## There is 1 tree of size {}".format(tree_size))
else:
  st.sidebar.write("## There are {} trees of size {}".format(len(prufer_codes), tree_size))
st.sidebar.write("## Of which, {} has Bounding Index {}".format(num_trees, bi))

st.sidebar.write("---")

num_samples = st.sidebar.number_input(
  "Number of Samples (up to {})".format(num_trees),
  min_value=min(5, num_trees),
  max_value=num_trees,
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

samples = random.choices(bi2prufer_codes[bi], k=num_samples)
for code in samples:
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
    for i in range(len(eccs_list)):
        if eccs_list[i]*len(eccs_list) < ecc_sum:
            bi = i
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
        "Prufer Code: [{}]\nAverage Eccentricity: {:.2f}".format(
            " ".join(map(str, code)),
            ecc_sum/tree_size
        )
    )
    st.pyplot(fig)
plt.clf()
    #dot = nx.nx_pydot.to_pydot(g)
    #st.graphviz_chart(dot.to_string
