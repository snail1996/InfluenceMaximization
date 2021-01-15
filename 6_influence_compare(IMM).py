import numpy as np
from igraph import *
import time
from tqdm import *#进度条显示
import matplotlib.pyplot as plt

def triggering(G, S, p, mc):          #igraph求影响力
    st = time.time()
    spread = []

    for _ in tqdm(range(mc)):

        new_nodes, RRS0 = S[:], S[:]

        while new_nodes:  # 求所选节点的激活节点集合

            B = G.neighborhood(vertices=new_nodes, order=1, mode="out")  # 新节点的邻居节点集合

            A = [item for subset in B for item in subset]  # 新节点的邻居节点集合，转化为一个列表

            success = np.random.uniform(0, 1, len(A)) < p

            temp = list(np.extract(success, A))  # 新节点以概率p随机激活邻居节点，每个节点只能激活一次

            RRS = list(set(RRS0 + temp))  # list(set())可以去除相同元素

            new_nodes = list(set(RRS) - set(RRS0))

            RRS0 = RRS[:]

        spread.append(len(RRS))

    return [np.mean(spread), time.time()-st]

S_CCA1 = [4400, 1458, 17915, 718, 6338, 15594, 638, 9796, 9514, 6218, 1450, 3433, 731, 9646, 1353, 40849, 2054, 18278, 6217, 5578]
S_DCA1 = [8312, 31011, 8700, 1458, 8764, 17915, 718, 17554, 20912, 12824, 53162, 22450, 13566, 5065, 22359, 9514, 12314, 6218, 1450, 5813]
S_CCA2 = [4400, 731, 18278, 3269, 18523, 28446, 14906, 36827, 41830, 116856, 37216, 57029, 73445, 22559, 34728, 83323, 121844, 110792, 67940, 47310]
S_DCA2 = [8312, 31011, 17915, 17554, 8802, 20128, 23741, 39631, 2243, 18577, 731, 132817, 17495, 28314, 11542, 13732, 108923, 4298, 8643, 47184]
g = Graph.Read_Edgelist("higgs-social_network_134574.txt", directed=True)
pr = 0.05
mc = 1000
k_l = range(1, 21)
i_CCA1 = []
i_CCA2 = []
i_DCA1 = []
i_DCA2 = []
#0.01概率的imm影响力
#i_IMM = [17.237, 32.951, 47.392, 62.892, 76.105, 89.528, 102.687, 113.926, 126.919, 135.773, 147.675, 158.039, 171.042, 181.31, 192.37, 199.978, 209.662, 220.297, 230.843, 241.303]
#0.05概率的imm影响力
i_IMM = [10011.414, 10065.082, 10089.89, 10126.783, 10159.098, 10148.714, 10164.229, 10201.763, 10219.956, 10247.786, 10266.077, 10283.634, 10321.188, 10324.258, 10328.144, 10351.961, 10387.034, 10408.184, 10410.526, 10438.712]

for k in k_l:
    print(k)
    # i_CCA1.append(IC(df, seed_sets1[:k], 0.05, 100)[0])
    # i_DCA1.append(IC(df, seed_sets2[:k], 0.05, 100)[0])
    t1 = triggering(g, S_CCA1[:k], pr, mc)
    i_CCA1.append(t1[0])
    t2 = triggering(g, S_DCA1[:k], pr, mc)
    i_DCA1.append(t2[0])
    t3 = triggering(g, S_CCA2[:k], pr, mc)
    i_CCA2.append(t3[0])
    t4 = triggering(g, S_DCA2[:k], pr, mc)
    i_DCA2.append(t4[0])

plt.figure()
plt.plot(k_l, i_CCA1, label='CCA1', marker='o', color='r', linestyle=':')
plt.plot(k_l, i_CCA2, label='CCA2', marker='+', color='y', linestyle='-')
plt.plot(k_l, i_DCA1, label='DCA1', marker='*', color='b', linestyle='--')
plt.plot(k_l, i_DCA2, label='DCA2', marker='D', color='k')
plt.plot(k_l, i_IMM, label='IMM', marker='^', color='g', linestyle='-.')
plt.legend()
plt.xlabel('seed_set size')
plt.ylabel('influence nodes')
plt.title('influence nodes vs seed_set size')
plt.show()