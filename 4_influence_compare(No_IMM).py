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

S_CCA1 = [4400, 1458, 17915, 718, 6338, 15594, 638, 9796, 9514, 6218, 1450, 3433, 731, 9646, 1353, 40849, 2054, 18278, 6217, 5578, 3430, 6212, 2025, 18523, 21752, 4306, 15592, 2465, 29383, 1105, 1946, 12839, 16911, 2470, 747, 18988, 2458, 9310, 28446, 14906, 13562, 5931, 1851, 2090, 9639, 10549, 1472, 5098, 63338, 16668]
S_DCA1 = [8312, 31011, 8700, 1458, 8764, 17915, 718, 17554, 20912, 12824, 53162, 22450, 13566, 5065, 22359, 9514, 12314, 6218, 1450, 5813, 20128, 2436, 39752, 6982, 20923, 19193, 39108, 23741, 7542, 3433, 2243, 23045, 5965, 29620, 438, 731, 39688, 9646, 56666, 22199, 4317, 938, 21229, 26122, 37642, 13612, 132817, 1353, 19678, 31016]
S_CCA2 = [4400, 731, 18278, 3269, 18523, 28446, 14906, 36827, 41830, 116856, 37216, 57029, 73445, 22559, 34728, 83323, 121844, 110792, 67940, 47310, 59703, 83477, 64126, 126654, 113044, 44378, 12229, 66978, 131810, 84488, 43074, 71087, 19119, 23505, 7393, 64014, 62587, 121180, 94873, 24745, 117742, 115751, 109029, 71486, 9726, 77736, 65619, 13716, 74585, 98522]
S_DCA2 = [8312, 31011, 17915, 17554, 8802, 20128, 23741, 39631, 2243, 18577, 731, 132817, 17495, 28314, 11542, 13732, 108923, 4298, 8643, 47184, 2251, 11927, 94542, 67153, 47735, 121844, 46282, 110792, 42619, 14906, 13716, 115751, 32811, 39028, 74585, 76982, 21909, 48197, 92221, 67003, 66978, 50884, 81676, 57898, 9390, 30101, 67940, 22363, 108352, 17128]

g = Graph.Read_Edgelist("higgs-social_network_134574.txt", directed=True)
pr = 0.05
mc = 1000
k_l = range(1, 51)
i_CCA1 = []
i_CCA2 = []
i_DCA1 = []
i_DCA2 = []
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
plt.plot(k_l, i_DCA2, label='DCA2', marker='^', color='g', linestyle='-.')
plt.legend()
plt.xlabel('seed_set size')
plt.ylabel('influence nodes')
plt.title('influence nodes vs seed_set size')
#plt.savefig('C:\\Users\\Snail\\Desktop\\picture\\test.png')
plt.show()