
from diversityStats.lib.uniprot_ec_dict import uniprot_ec_dict
from diversityStats.lib.gini_simpson import gini_simpson_dict, gini_simpson_value
from decimal import Decimal
from scipy.cluster.hierarchy import fcluster, linkage
from collections import defaultdict
import numpy as np
import json
from collections import defaultdict
import matplotlib.pyplot as plt
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
import seaborn as sns
from copy import deepcopy


_MAT = "./PF00171_mat.npy"
_HEAD = "./PF00171_headings.json"
_ANN = "./PF00171_annotation.tab"


def sep_indices(val):

    set_indices = []
    non_indices = []

    for i in range(0, len(val)):
        if val[i] == 1:
            set_indices += [i]
        else:
            non_indices += [i]

    return (set_indices, non_indices)


def score_sol(mat, sol):

    score = 0
    (set_indices, non_indices) = sep_indices(sol)

    sim_list = []

    for i in range(0, len(set_indices)):
        for j in range(i + 1, len(set_indices)):
            sim = ((1 - mat[set_indices[i], set_indices[j]])) - (mat[set_indices[i], set_indices[j]]/(1 - mat[set_indices[i], set_indices[j]]))
            score += sim

            sim_list.append(mat[set_indices[i], set_indices[j]])

    return score, sim_list


def get_ec_subset(subset, ac_to_ec):

    ec_dict = defaultdict(int)

    for ac in subset:
        ec_str = ac_to_ec[ac]
        ec_lst = ec_str.split(';')

        for ec in ec_lst:
            ec = ec.strip()
            ec_dict[ec] += 1

    return ec_dict


def initialise_matrix(dist_file):

    dist = np.load(dist_file)

    for i in range(0, len(dist)):
        dist[i][i] = np.nan
        for j in range(i+1, len(dist)):
            if dist[i][j] == 1:
                dist[i][j] = 0.99
                dist[j][i] = 0.99

    return np.asmatrix(dist)


def initialise_headings(heading_file):

    with open(heading_file) as heading:
        headings = json.loads(heading.read())

    temp_dict = {}

    for i in range(0, len(headings)):
        temp_dict[i] = headings[i]

    return temp_dict


def random_solution(length, num_picked):

    arr = np.array([0] * (length - num_picked) + [1] * num_picked)
    np.random.shuffle(arr)
    return list(arr)


def get_n(Z, n, mat):

    tmax = Decimal("1000")
    tmin = Decimal("0")
    tcurrent = Decimal("500")

    clustering = fcluster(Z, tcurrent)
    result = len(set(clustering))

    while result != n:
        if result < n:
            tmax = tcurrent
            tcurrent = (tmin + tcurrent) / 2
        else:
            tmin = tcurrent
            tcurrent = (tmax + tcurrent) / 2

        clustering = fcluster(Z, tcurrent)
        result = len(set(clustering))

        print(f"{tcurrent}: {result}")

    clust_dict = defaultdict(list)

    for i in range(len(clustering)):
        curr_clust = clustering[i]
        clust_dict[curr_clust].append(i)

    sol = np.zeros(len(mat))

    for indices in clust_dict.values():
        chosen_ind = np.random.choice(indices)
        sol[chosen_ind] = 1

    return sol


def main():

    mat = initialise_matrix(_MAT)
    head = initialise_headings(_HEAD)
    ac_to_ec, ec_to_ac = uniprot_ec_dict(_ANN, 2)

    # pf00171
    custom_subset = ['Q6F9F7', 'A5GSH0', 'P0DPF0', 'B8I6T0', 'Q2FWX9', 'P43503', 'Q8GAI8', 'P46329', 'Q55167', 'P55653', 'Q9A777', 'H2IFE7', 'O05619', 'O32507', 'B2V9F3', 'Q5E2G9', 'Q985M6', 'P86808', 'B8DCT8', 'Q6F9G0', 'Q738L2', 'Q9I6C8', 'B7MVM5', 'Q2SKP1', 'P80668', 'Q3B2U3', 'Q55585', 'B3QPW0', 'P25553', 'P28810', 'Q2BN77', 'Q92YD2', 'Q8GAK7', 'P39616', 'Q5X4K4', 'Q2YV11', 'A5YBJ3', 'Q7VBM1', 'A6TUA0', 'Q59702', 'P42269', 'P23883', 'P23105', 'C0R0B8', 'B2VJX8', 'Q6G4Z0', 'A1WGI4', 'H8ZPX2', 'Q84DC3', 'P94358', 'Q2G9T9', 'Q7VI05', 'P0C1E0', 'Q3IC91', 'P76149', 'Q8CNI5', 'H1ZV37', 'P38947', 'B0S9A5', 'P12693', 'B2IZ89', 'Q47UQ0', 'E1V7V8', 'Q1JUP4', 'Q4JWT3', 'Q6LTX2', 'A1VYR7', 'A5CXP4', 'Q0K845', 'O86447', 'Q1GV29', 'A8EVN0', 'Q7N2G9', 'A5FYS4', 'A2C148', 'Q9I702', 'Q79EM7', 'Q9AHG1', 'P33008', 'Q165Y8', 'B2S2U7', 'O69497', 'Q9KWS5', 'B1XDF5', 'Q03ZF1', 'P0A391', 'Q5WH11', 'Q59931', 'B0U208', 'A2BQ71', 'B0T8I8', 'Q8CUQ4', 'Q72NQ9', 'Q4QJW6', 'A3DC22', 'B4U8A0', 'Q92UV7', 'Q04FB2', 'Q12SX8', 'Q4L919']
    #pf00155
    # custom_subset = ['Q8DTM1', 'Q058A6', 'C6C2Z3', 'Q6LT75', 'B1N009', 'P9WQ88', 'P39643', 'Q9I468', 'B8ERL9', 'P43089', 'P72173', 'C3KVX5', 'P95468', 'Q6MDE0', 'Q253K9', 'B6YRL2', 'P9WPZ4', 'B9KDN6', 'P43336', 'P04693', 'Q8DM42', 'Q08432', 'B0TY45', 'Q5SHW0', 'Q8XZC3', 'Q68XV9', 'Q9CJU0', 'Q92G23', 'Q9KM65', 'Q9Z856', 'Q8ABA8', 'Q7UZZ3', 'Q81V80', 'B1WY56', 'Q5F6R6', 'A6L8U2', 'P00509', 'B4RFX5', 'Q8KZM9', 'P0A959', 'B7VH15', 'A0L3L7', 'Q4FP52', 'Q3S8P9', 'P28735', 'Q02636', 'C0QFJ4', 'Q8R5U4', 'P96847', 'O25320', 'P63503', 'Q5ZW88', 'P77806', 'P74770', 'B8D707', 'P97084', 'Q72LL6', 'Q9RRM7', 'A4SPR6', 'Q6D3C0', 'Q84CG1', 'B4U9L1', 'Q7VL09', 'P09053', 'A5N7Q7', 'Q93QC6', 'Q7W9I4', 'P44425', 'Q492K2', 'B2A250', 'Q3ARM7', 'A7N6R9', 'Q55128', 'O07587', 'P36692', 'Q3BYN0', 'Q9L6I2', 'Q89AX7', 'B5YFU5', 'Q06965', 'A8MEH2', 'P37419', 'Q9CBM9', 'Q5FRR4', 'Q89AK6', 'Q72PG3', 'P16524', 'P21633', 'P63499', 'A5INE2', 'O87320', 'B8CX89', 'O52815', 'Q0BBD6', 'P0A4X5', 'A0PXP5', 'Q02135', 'A5CVR5', 'P36570', 'B8J3V0']

    #pf04055
    # custom_subset = ['B3E599', 'A5G2D2', 'Q057Q1', 'B8FS78', 'Q9K864', 'Q7NCE3', 'Q8CJT5', 'Q185C5', 'Q97L63', 'Q7NIT2', 'Q2SWB9', 'A9FD89', 'P71011', 'A4XGB8', 'Q82K95', 'A1W1T3', 'O25376', 'Q057G5', 'Q9EYN8', 'B5QX73', 'B7GQG0', 'P39409', 'Q4FNN5', 'Q38HX2', 'Q057Q7', 'A0A1C7D1B7', 'B0CDZ6', 'P43751', 'P09825', 'P32131', 'Q0RD46', 'Q55373', 'P10390', 'Q8K9D9', 'B2RH08', 'P74132', 'Q8KBK9', 'A0Q2E1', 'A6H1N2', 'B8ENI9', 'A7ZE07', 'O83293', 'O67826', 'B0VQD7', 'Q8KC85', 'Q8YR77', 'B2GLQ7', 'B9M4F4', 'Q2RSY6', 'Q49573', 'Q89ZC3', 'P17434', 'A9EPV3', 'P20714', 'A8Z642', 'Q17XY7', 'Q8KCU0', 'Q9K0Q5', 'O34162', 'Q8KFK8', 'B1Y6D6', 'Q3J561', 'Q8EUX4', 'P9WJ78', 'B1IL14', 'Q81G67', 'Q02550', 'P51008', 'P73667', 'A0A384LP51', 'Q55914', 'Q53U14', 'P45097', 'Q8D1Y5', 'A0LV48', 'O33506', 'Q2JRI4', 'P55477', 'Q72DS4', 'A0RIB6', 'Q0TTH1', 'Q44634', 'A8Z609', 'B1ZVM5', 'P69848', 'C9XIS7', 'Q1IHK7', 'P75794', 'A9A0B5', 'Q6MED6', 'A6LSR6', 'P24427', 'Q8RHX4', 'Q30XT6', 'A0A069AMK2', 'O87941', 'Q9S498', 'Q8DII8', 'Q1GV98', 'A9CF16']

    # pf00155 CSN
    # custom_subset = ['Q8XZC3', 'Q11VM5', 'B7GHW7', 'Q56232', 'Q3AMU5', 'Q8EFB2', 'O52815', 'C0QFJ4', 'Q72LL6', 'Q84CG1', 'P09053', 'Q04UL5', 'Q39M27', 'Q02635', 'Q8Z8H8', 'Q2JPM4', 'Q9RI00', 'B9M384', 'Q8DM42', 'Q7NL03', 'P36692', 'B8DJJ6', 'Q9L6I2', 'Q9I468', 'P28735', 'Q93QC6', 'P63501', 'A0PVN0', 'Q30TC9', 'P61001', 'Q311Z4', 'B2I5Y0', 'Q8KZM9', 'Q7N6Q6', 'Q7M7Y6', 'P26505', 'B7KD70', 'B9MDV4', 'B8ZR84', 'A0RIB9', 'P9WML7', 'Q9JTH8', 'Q87QN5', 'Q9JYH7', 'Q8DJ97', 'Q5F7D7', 'B3ECG2', 'Q058A6', 'A3MNG3', 'Q5H0L0', 'Q8Z5J9', 'O07587', 'Q08432', 'Q89AX7', 'Q1GP30', 'Q9HVX0', 'Q492K2', 'P61004', 'Q1B1Z8', 'P73807', 'Q02636', 'P74861', 'Q9K0U0', 'Q9CBM9', 'P21633', 'Q3S8P9', 'Q3SLX9', 'P44425', 'P43336', 'Q10VS0', 'Q12F38', 'Q119K2', 'Q7P0F4', 'B1LYP9', 'Q12D74', 'B8HW95', 'Q0BUV6', 'Q7VL09', 'P0A959', 'Q1QYD6', 'Q1LS75', 'Q68VS3', 'Q89AK6', 'P36570', 'Q824A4', 'Q9Z856', 'Q3A3Z8', 'P18079', 'Q9KM65', 'Q04QW8', 'B0B7W0', 'Q8G4S8', 'Q8D8N0', 'Q9ZLN3', 'Q48ED0', 'B9K9R9', 'A9HVE9', 'Q5Z3C0', 'Q9KSX2', 'A9R3C9']

    # pf00155 kmedoids
    # custom_subset = ['A1VUJ6', 'Q02636', 'A8ZUS7', 'A7N6R9', 'B0JQZ0', 'B4RFX5', 'Q9L6I2', 'Q9KM65', 'B8J3V0', 'Q9CBM9', 'Q9KSZ3', 'P22806', 'O07587', 'B0U6J0', 'Q9ZHE5', 'P00509', 'Q5F6R6', 'P09053', 'Q7VL09', 'Q84I52', 'P0A961', 'A1K6Q1', 'Q7WDY3', 'Q02YW3', 'P36692', 'Q47829', 'A0PP02', 'P58661', 'Q3SLX9', 'Q2NUJ6', 'Q87QN5', 'Q6D3C0', 'A8AJ11', 'P37419', 'Q56114', 'P95468', 'A6L2V8', 'Q3A3Z8', 'B2VBT8', 'B3PI88', 'P9WQ91', 'P63499', 'Q0A5W2', 'Q08432', 'Q7VQW9', 'B2SWS7', 'Q7W9I4', 'Q5L6M0', 'Q9K625', 'Q62GE0', 'A2SD53', 'B5EGX2', 'B8CX89', 'B2IH50', 'C0QFJ4', 'B7GHW7', 'Q73KM3', 'Q3Z8H5', 'Q4K4T3', 'B2FLM5', 'O67857', 'B4ESU4', 'Q89AK6', 'P57202', 'Q81I05', 'Q68VS3', 'B5ELF7', 'Q55128', 'Q06965', 'A1AQT1', 'A1ITK1', 'Q119K2', 'Q5SHW0', 'Q3BYN0', 'Q9ZLN3', 'B2A250', 'A7ZY32', 'A0QHJ9', 'Q058A6', 'Q5KY23', 'A9VG56', 'B5BC30', 'P44422', 'B5YFU5', 'P97084', 'Q5NZF5', 'P26505', 'Q1LS75', 'B8ERL9', 'A7HP29', 'Q6LPR3', 'P18079', 'Q2JS04', 'B8F713', 'Q1LT68', 'A1WVM6', 'A9HJ57', 'A7GSE1', 'Q8DJ97', 'Q39Z65']

    # pf04055 kmedoids

    # custom_subset = ['B5ECN1', 'Q057G5', 'Q8KFK8', 'P45097', 'Q02550', 'Q8K9D9', 'Q9K0Q5', 'Q0TTH1', 'A6H1N2', 'P36569', 'P73667', 'Q609V2', 'P55477', 'A0Q2E1', 'Q8XMQ3', 'Q7NCE3', 'A0L3M0', 'E3PRJ8', 'P9WJ78', 'A8GBC5', 'Q8KD71', 'Q72DS4', 'Q55373', 'A8Z642', 'Q7NIT2', 'Q2RSY6', 'B7ID25', 'P57373', 'A8Z609', 'A7NPY6', 'Q9ZLH0', 'Q9X758', 'A0LV48', 'B1ZQZ5', 'Q44634', 'P74132', 'P72811', 'A7GSD8', 'P30140', 'A0L509', 'P32131', 'Q18CP3', 'Q81G67', 'A8G2A6', 'A9NEU7', 'Q72DM4', 'B8DPP9', 'Q6F9I9', 'Q55914', 'Q56211', 'C1C6B0', 'Q057Q1', 'A9CF16', 'B4ESU3', 'Q49573', 'Q3JCN4', 'B6IQ36', 'Q0RDQ8', 'A0A384LP51', 'Q53U14', 'Q8DII8', 'Q51676', 'Q38HX2', 'Q8EVK0', 'Q2RT68', 'Q7U3H2', 'B1HQE6', 'Q899M1', 'B3DYX1', 'Q8YR77', 'Q057Q7', 'Q7VFE6', 'Q3ASD4', 'B3E599', 'Q2JRI4', 'B6JNP9', 'P75794', 'A7I414', 'Q3AVP9', 'Q8DLC2', 'A8L6D8', 'Q2J5A7', 'A0PKZ7', 'Q3AGB7', 'B8ENI9', 'Q7MT97', 'Q7U7Q2', 'B4SHB4', 'Q89AK8', 'Q8D1Y5', 'B3CQP5', 'Q74CF3', 'Q0BSW6', 'P57357', 'A6LUQ4', 'B1JE55', 'A5I5U3', 'A9A0B5', 'B6JNS0', 'A6LD84']

    # pf00171 kmedoids
    # custom_subset = ['Q8GAI8', 'Q6HJ19', 'P46329', 'Q2BN77', 'Q1QTQ7', 'Q5ZUT5', 'Q65IX1', 'Q5X4K4', 'Q63BL0', 'Q5WVZ4', 'Q2FWX9', 'Q9I702', 'Q3IC91', 'Q738L2', 'A0RDW1', 'B2IZ89', 'A9VK31', 'Q0I8Z0', 'A7GPH3', 'Q0K845', 'P42269', 'A2BVQ3', 'Q2SKP1', 'C1EZ15', 'B5F7J0', 'Q9A777', 'A8GHZ8', 'P9WNX8', 'Q59931', 'Q59702', 'A9VF06', 'A8AHD8', 'A8G3V6', 'Q3C1A6', 'Q81QR5', 'A2C148', 'Q84DC3', 'P55653', 'A5TYV9', 'Q7VBM1', 'P94358', 'P86808', 'Q55585', 'P39616', 'B0CFL0', 'P38947', 'B9IZZ7', 'P33008', 'A1V675', 'B1KDF2', 'P53000', 'P25553', 'C1CRW1', 'O69497', 'P80668', 'O86447', 'B7HR31', 'Q4L919', 'A8EVN0', 'Q4L803', 'B0S9A5', 'P0DPF0', 'P0A390', 'P0A391', 'B2V9F3', 'B7H597', 'Q7M8Z4', 'B0KR47', 'Q9I6M5', 'C0R0B8', 'Q112S1', 'B2VJX8', 'Q7A1Y7', 'Q7A825', 'Q6GKD8', 'A9EN94', 'Q839W3', 'B9L9A2', 'A1U5W8', 'Q6GCV9', 'Q8ES27', 'O05619', 'O06478', 'A2BQ71', 'A2CAS7', 'A7ZAI1', 'Q30SG0', 'Q1GV29', 'A6WST7', 'Q5WH11', 'Q99X54', 'B7IW48', 'Q8CNI5', 'E1V7V8', 'Q8Y9Y4', 'Q7NXX7', 'Q5L025', 'Q8EMV4', 'Q723T1', 'B8DCT8']

    # pf00155 kmedoids size 200

    # custom_subset = ['P37419', 'P63499', 'A9VG56', 'B4ESU4', 'Q9Z856', 'A8ZUS7', 'Q9L6I2', 'Q9KM65', 'Q9HUI9', 'Q9CBM9', 'Q93QC6', 'B1Y500', 'A7HP29', 'Q8Z8H8', 'B2VBT8', 'Q8KDS8', 'B0UKC8', 'Q8DTM1', 'B8ERL9', 'Q84CG1', 'Q82WA8', 'Q824A4', 'Q7VA14', 'A7Z5B4', 'B0U6J0', 'Q72LL6', 'Q6MDE0', 'B7ULX3', 'B2SWS7', 'Q68VS3', 'B3PI88', 'Q5SHW0', 'A2SD53', 'Q5L6M0', 'Q56114', 'Q55128', 'Q9K625', 'Q4UJV4', 'Q6D3C0', 'A7FKM8', 'Q7VQW9', 'Q02YW3', 'Q2JS04', 'Q2JLL9', 'Q253K9', 'Q7W9I4', 'Q1RIV2', 'Q2W3L2', 'B2IH50', 'P22806', 'A8AJ11', 'Q08432', 'Q06191', 'A7ZY32', 'P57202', 'A0QHJ9', 'Q04YV8', 'Q04UL5', 'Q02636', 'P9WQ91', 'P9WQ90', 'A1ITK1', 'Q2Y9Y8', 'B7K0L9', 'P55683', 'P97084', 'P96847', 'P95468', 'P44422', 'P77806', 'Q73KM3', 'P74861', 'P72173', 'A6L2V8', 'B2UBJ3', 'B1XL23', 'P58661', 'B7NA75', 'Q0A5W2', 'P44425', 'P43336', 'B5YFU5', 'P36692', 'P28735', 'P23034', 'P16524', 'P0A961', 'P0A960', 'P0A959', 'P09053', 'P04693', 'P00509', 'Q72PG3', 'O85746', 'Q3JMZ7', 'Q84I52', 'O52815', 'Q3A3Z8', 'O07587', 'B2FLM5', 'C6C2Z3', 'A5INE2', 'B8CX89', 'Q81I05', 'B6YRL2', 'B5EGX2', 'B2A250', 'Q2NUJ6', 'A9HVE9', 'Q3BYN0', 'B0SEH8', 'B0CDH5', 'Q21FY4', 'B5F073', 'Q87QN5', 'A8ZXV5', 'A7N6R9', 'A9HJ57', 'Q9ZLN3', 'C4Z4Y1', 'A1ATI6', 'Q74GT3', 'Q8PQD8', 'B3E933', 'B5ELF7', 'O25320', 'P9WPZ4', 'A8MEH2', 'Q02635', 'A1VUJ6', 'P21633', 'Q9PK04', 'B6ESC6', 'Q82DR2', 'P43089', 'C4ZG66', 'O86459', 'Q8D8N0', 'B3Q8Z5', 'Q8DJ97', 'A3PMF8', 'A7GSE1', 'B7VH15', 'B4U9L1', 'P18079', 'Q06965', 'Q058A6', 'Q9I468', 'A2CC97', 'Q04512', 'B8HJY4', 'A5UCE4', 'Q3MDN5', 'A1K6Q1', 'Q7VL09', 'A8MEX7', 'Q7NDX4', 'Q795M6', 'A1TKZ0', 'A5GIN1', 'Q72BI1', 'Q3SLX9', 'Q39XE0', 'Q8YP73', 'Q1LT68', 'A5FZN8', 'Q60013', 'Q5KY23', 'Q5NZF5', 'Q9KSZ3', 'Q59228', 'A9BV10', 'Q55828', 'B8J3V0', 'Q89AK6', 'A1AQT1', 'O67781', 'Q8DH57', 'Q3AW44', 'Q47829', 'Q4K4T3', 'O31665', 'Q2RK33', 'Q119K2', 'A1WVM6', 'Q6LPR3', 'Q3ZXC8', 'Q9ZHE5', 'Q3AC10', 'B8F713', 'Q9K0U0', 'A1KVF6', 'P26505', 'Q5F6R6', 'B4RFX5', 'A5N7Q7', 'C0R1Z0', 'B4RP93', 'P08080', 'Q6GJB8']

    # pf00171 kmedoids size 200
    # custom_subset = ['Q6HJ19', 'B8DCT8', 'P46329', 'Q2BN77', 'Q65D00', 'P43503', 'Q65IX1', 'Q9KWS5', 'Q63BL0', 'A0PN13', 'Q9KAH5', 'Q9I702', 'A0QMB9', 'Q9I6M5', 'A0RDW1', 'A0REB5', 'P42412', 'Q9I6C8', 'A9VMS6', 'P42329', 'P42269', 'A8EVN0', 'A9EN94', 'Q9A7W2', 'A1KF54', 'Q9A777', 'Q99X54', 'Q1QTQ7', 'Q99SD6', 'B5FBM2', 'A9VF06', 'Q92UV7', 'Q92EQ7', 'Q3C1A6', 'Q63B74', 'B0S9A5', 'P42236', 'P55653', 'A9MQY3', 'P71016', 'Q6D6Y7', 'B7VLI4', 'Q2FF06', 'P39616', 'Q6F9F7', 'P38947', 'B4ECZ8', 'P33008', 'P28810', 'Q8Y9Y4', 'Q488Y0', 'P25553', 'B6EML5', 'P25526', 'A2BVQ3', 'P23883', 'Q839W3', 'P19059', 'P12693', 'Q8NVG4', 'Q6F9G0', 'P0DPF0', 'B1JCH3', 'P0DOV9', 'P76149', 'Q49Z69', 'Q8GAK7', 'Q8GAI8', 'B2VJX8', 'B7H597', 'A4IPB2', 'A4IPF5', 'P0A391', 'P0A390', 'Q4L803', 'Q4L919', 'B7HR31', 'O86447', 'Q4V1F6', 'A4SK35', 'Q8ES27', 'Q8EMV4', 'P80668', 'Q6FCQ0', 'O69763', 'O69497', 'Q8DCS9', 'Q2YV11', 'Q7V8C3', 'Q2YUN1', 'A4WAR9', 'B7IW48', 'Q8CNI5', 'Q8CN24', 'O32507', 'Q7NXX7', 'Q88RC0', 'Q1GV29', 'B9IZZ7', 'O06478', 'O05619', 'H2IFE7', 'H1ZV37', 'Q15Y60', 'E1V7V8', 'Q55585', 'P86808', 'Q3IC91', 'C6DD82', 'P94358', 'A5ICK3', 'Q2SKP1', 'P94428', 'A5TYV9', 'Q84DC3', 'Q141D3', 'Q81QR5', 'Q81QB6', 'Q81DR5', 'A8AHD8', 'Q6G7I8', 'Q59702', 'Q59931', 'Q5E2G9', 'P9WNX8', 'P9WNX9', 'A5YBJ3', 'Q5HE78', 'C3LRS7', 'B0T8I8', 'C1CRW1', 'Q7U2I0', 'Q5HJK3', 'C1KZ99', 'Q5HLA3', 'Q5HMA0', 'Q7N2G9', 'Q5KYK0', 'Q6GCV9', 'H8ZPX2', 'Q5KYR4', 'Q0K845', 'A6WST7', 'Q5L025', 'B1KDF2', 'A7BJC4', 'Q2G9T9', 'Q6GEV3', 'A7GPH3', 'Q6GKD8', 'Q7A825', 'Q7A4D8', 'Q7A1Y7', 'Q2FK94', 'A7MX96', 'Q5QVX3', 'Q79EM7', 'Q73TP5', 'A7ZAI1', 'Q738L2', 'Q2G1J0', 'Q2FWX9', 'Q5WH11', 'Q2FWD6', 'Q6HIK3', 'Q5WKZ1', 'Q5WVZ4', 'A8GHZ8', 'Q5X4K4', 'Q723T1', 'Q5ZUT5', 'C0R0B8', 'A8FDV4', 'A8G3V6', 'B2IZ89', 'P53000', 'A2C148', 'Q2SXN9', 'Q30SG0', 'A0KN18', 'Q3IFT7', 'Q5PHV8', 'Q0I8Z0', 'B0KR47', 'Q7MH21', 'Q7M8Z4', 'Q62LN5', 'Q7V293', 'A5F529', 'Q8Z747', 'A2CAS7', 'Q9KNW4', 'A1U5W8', 'Q87L22', 'Q3K885', 'B5F7J0', 'B1JYT5', 'Q8ZPC9', 'Q31BU4', 'A9MYQ4']

    #  pf04055 kmedoids size 200
    # custom_subset = ['Q56211', 'Q51676', 'A0A384LP51', 'A8Z642', 'Q72DS4', 'Q2NYZ2', 'B7ID25', 'Q9ZLH0', 'A6H1N2', 'A0A1C7D1B7', 'Q38HX2', 'B8D983', 'A0RQM9', 'Q74CF3', 'A1K1U7', 'Q30Y73', 'Q0BSW6', 'Q89AK8', 'P07748', 'A5N854', 'Q8D2A1', 'Q8KD71', 'Q609V2', 'P09825', 'Q5FFY2', 'Q2RSY6', 'C1DEW5', 'Q30XT6', 'Q0C0U1', 'B6JNS0', 'Q7NCE3', 'Q9X758', 'A0LV48', 'A5FY82', 'Q44634', 'Q728P5', 'A5CCM4', 'P0A1E2', 'Q7NFJ9', 'Q97I18', 'A7I414', 'A0L509', 'Q8DLC2', 'B1IL14', 'B2FUS8', 'A3PJL9', 'Q3AGB7', 'Q6F9I9', 'A6LD84', 'A9EPV3', 'A0Q2E1', 'Q2J5A7', 'A4IM49', 'O31677', 'Q49573', 'A5EVN1', 'B3EGT4', 'A8L6D8', 'B9J9I5', 'Q12F39', 'O33506', 'P57373', 'P20627', 'A5GW06', 'O34162', 'A9CF16', 'Q81G67', 'B5ECN1', 'A1WE19', 'Q18CP3', 'Q1QRH1', 'P32131', 'Q92CY2', 'Q8K9P5', 'B1X7B1', 'A1AWE7', 'Q84F14', 'B9KEA4', 'B1XXL6', 'A9FD62', 'B3E599', 'Q4ZMC1', 'B4SHB4', 'Q5SGZ3', 'B8J0C6', 'Q0TTH1', 'Q7MT97', 'B6YRD1', 'Q2GG00', 'B1IHH9', 'Q057Q1', 'Q9K0Q5', 'Q3JCN4', 'Q53U14', 'Q8K9D9', 'Q55373', 'B0SGV6', 'A9BGV7', 'Q89AK5', 'A9KG71', 'Q55914', 'Q02550', 'Q01060', 'Q8KFK8', 'B6JNP9', 'P9WJ78', 'B7VAB6', 'P45097', 'P95651', 'Q3YRG6', 'Q8D1Y5', 'A1B1Z0', 'B2RMI0', 'P75794', 'P74132', 'P73667', 'B5EMR4', 'Q5HJF3', 'Q057G5', 'Q97L63', 'P55477', 'Q2RT68', 'P71517', 'Q899M1', 'Q8XMQ3', 'Q97D55', 'C1C6B0', 'A9FD89', 'Q7NCF5', 'A9WFY6', 'O67886', 'P59749', 'B1JE55', 'B2A3C0', 'Q7NIT2', 'Q0RDQ8', 'B2SSW3', 'Q8KBK9', 'Q8KCU0', 'A7NPY6', 'B9L7E1', 'O87941', 'Q7U3H2', 'A2C8D5', 'E3PRJ8', 'Q8YR77', 'P30140', 'Q3IGS6', 'A9A0B5', 'A5VUG1', 'B2J1M7', 'Q9L3B0', 'Q0STS9', 'Q7NDA8', 'B8EIR0', 'Q3ASD4', 'B1ZQZ5', 'Q2IUT3', 'A1ULP7', 'Q2NAD7', 'B2TL73', 'Q5LTB4', 'Q7U7Q2', 'Q9S498', 'Q8CJT5', 'A7I1A0', 'B4ESU3', 'B3QIL2', 'B1Y6D6', 'Q895H1', 'B1HQE6', 'Q8EVK0', 'P72811', 'A8GBC5', 'B3ES11', 'Q2JRI4', 'Q5HSB2', 'Q3AVP9', 'P36569', 'Q31PU7', 'Q3YRT2', 'A9KK15', 'B3DYX1', 'Q057Q7', 'B6IQ36', 'Q72DM4', 'A6LUQ4', 'Q8DII8', 'A2SLX7', 'A7H0I4', 'A8G2A6', 'Q5P7B0', 'A1KJ05', 'B8DPP9', 'B4S808', 'Q5YZ59', 'Q7VFE6', 'Q39VB0', 'P57357', 'A5G3Q7']

    # pf04055 kmedoids size 50
    # custom_subset = ['Q2IUT3', 'Q057Q7', 'Q8KFK8', 'P45097', 'Q02550', 'Q8K9D9', 'Q7VFE6', 'E3PRJ8', 'A6H1N2', 'Q057G5', 'A6LD84', 'Q7MT97', 'P55477', 'A0Q2E1', 'Q8XMQ3', 'Q7NCE3', 'Q53U14', 'B6JNP9', 'A7GH77', 'A0LV48', 'Q72DM4', 'A9CF16', 'A3PJL9', 'Q55373', 'Q81G67', 'Q18CP3', 'Q7U7Q2', 'P32131', 'Q2JRI4', 'B1ZVM5', 'Q9ZLH0', 'Q9X758', 'Q0RDQ8', 'P73667', 'A8Z642', 'Q8D1Y5', 'B4SHB4', 'P74132', 'P30140', 'Q3AGB7', 'Q2RSY6', 'Q44634', 'Q057Q1', 'Q9K864', 'A6LUQ4', 'A7NPY6', 'A0A1C7D1B7', 'B3E599', 'Q8YR77', 'B1ZQZ5']

    # pf00155 kmedoids size 50
    # custom_subset = ['A2SD53', 'A5CVR5', 'A8ZUS7', 'A7N6R9', 'Q5NZF5', 'Q9ZLN3', 'Q9L6I2', 'Q9KM65', 'Q68VS3', 'Q9CBM9', 'P18079', 'Q7W9I4', 'Q5F6R6', 'P22806', 'B5F073', 'A7GSE1', 'Q3BYN0', 'P09053', 'Q2NUJ6', 'Q3A3Z8', 'Q058A6', 'Q119K2', 'P28735', 'Q06965', 'P36692', 'Q89AK6', 'C3KVX5', 'B2IH50', 'C6E9Q7', 'B4ESU4', 'Q84I52', 'B5YFU5', 'Q6LPR3', 'Q2JS04', 'Q3SLX9', 'P95468', 'B8ERL9', 'B2A250', 'Q73KM3', 'A1KVF6', 'B3PI88', 'P57202', 'Q0A5W2', 'Q7VL09', 'B2FLM5', 'B2VBT8', 'B8CX89', 'A7HP29', 'B8F713', 'Q6AL81']

    # pf00171 kmedoids size 50
    # custom_subset = ['Q7M8Z4', 'A1U5W8', 'A2BQ71', 'B9L9A2', 'A9VK31', 'Q5ZUT5', 'Q112S1', 'Q5X4K4', 'Q7NXX7', 'Q5WVZ4', 'Q2FWX9', 'P0DPF0', 'Q9I6M5', 'A8EVN0', 'Q4L919', 'Q839W3', 'C1CRW1', 'A1V675', 'A7GPH3', 'A2BVQ3', 'P42269', 'A8AHD8', 'C0R0B8', 'B2VJX8', 'Q2SKP1', 'A6WST7', 'P53000', 'B1KDF2', 'Q30SG0', 'A2CAS7', 'P33008', 'B0S9A5', 'A8G3V6', 'A9EN94', 'B0CFL0', 'P39616', 'Q84DC3', 'B0KR47', 'P86808', 'Q0I8Z0', 'Q3C1A6', 'C1EZ15', 'Q8GAI8', 'A2C148', 'Q1QTQ7', 'P80668', 'Q3IC91', 'Q5LRY6', 'B2IZ89', 'P9WNX8']

    # pf00171 kmedoids size 150
    # custom_subset = ['Q8GAI8', 'Q6HJ19', 'P46329', 'Q2BN77', 'Q65D00', 'Q5ZUT5', 'Q65IX1', 'Q5X4K4', 'Q63BL0', 'Q5WVZ4', 'Q2FWX9', 'Q9I702', 'A0QMB9', 'Q738L2', 'A0RDW1', 'Q73TP5', 'Q79EM7', 'Q9I6C8', 'A7GPH3', 'Q0K845', 'P42269', 'A6WST7', 'A1V675', 'Q7U2I0', 'A1KF54', 'Q9A777', 'P9WNX9', 'P9WNX8', 'Q59931', 'Q59702', 'A9VF06', 'Q92UV7', 'A8G3V6', 'Q3C1A6', 'Q81QR5', 'B0S9A5', 'Q84DC3', 'P55653', 'A5TYV9', 'B2IZ89', 'P94358', 'P86808', 'Q55585', 'P39616', 'H2IFE7', 'P38947', 'B9IZZ7', 'P33008', 'B2VJX8', 'B1KDF2', 'Q0I8Z0', 'P25553', 'Q7M8Z4', 'O69497', 'P80668', 'O86447', 'B7HR31', 'Q4L919', 'A8EVN0', 'Q4L803', 'A2BVQ3', 'P0DPF0', 'P0A390', 'P0A391', 'P76149', 'B7H597', 'Q8GAK7', 'Q2YV11', 'E1V7V8', 'A7ZAI1', 'Q112S1', 'Q2FK94', 'Q7A1Y7', 'Q7A825', 'Q6GKD8', 'B5F7J0', 'A8AHD8', 'Q2G1J0', 'Q5L025', 'Q6GCV9', 'Q8ES27', 'Q8EMV4', 'A8FDV4', 'Q5HMA0', 'B1JCH3', 'C1KZ99', 'P43503', 'B8DCT8', 'Q5HJK3', 'Q5WH11', 'Q723T1', 'B7IW48', 'Q8CNI5', 'P53000', 'Q99X54', 'Q8Y9Y4', 'Q92EQ7', 'Q1GV29', 'Q9KWS5', 'O06478', 'O05619', 'Q1QTQ7', 'Q839W3', 'C1CRW1', 'Q5WKZ1', 'Q6F9F7', 'Q9I6M5', 'Q7A4D8', 'A9VMS6', 'Q2G9T9', 'Q9KAH5', 'A9VK31', 'P94428', 'Q2FWD6', 'Q5HE78', 'Q99SD6', 'Q49Z69', 'Q88RC0', 'Q2FF06', 'P42329', 'Q81QB6', 'Q63B74', 'Q6HIK3', 'A0REB5', 'Q5QVX3', 'Q30SG0', 'Q3IC91', 'P42236', 'Q8NVG4', 'Q6GEV3', 'Q2YUN1', 'P19059', 'A2CAS7', 'P25526', 'A9EN94', 'A2C148', 'C0R0B8', 'Q5KYK0', 'Q6G7I8', 'Q2SKP1', 'Q4V1F6', 'Q9A7W2', 'Q5KYR4', 'B9L9A2', 'A4IPF5', 'B6EML5', 'Q6D6Y7', 'A1U5W8', 'A2BQ71', 'A4IPB2']

    # pf00155 kmedoids size 150
    # custom_subset = ['A7HP29', 'P63499', 'B5ELF7', 'Q9ZCB8', 'Q21FY4', 'A8ZUS7', 'Q9L6I2', 'Q9KM65', 'Q9HUI9', 'Q9CBM9', 'Q93QC6', 'Q0A5W2', 'A1ITK1', 'A7ZY32', 'Q3A3Z8', 'Q8KDS8', 'Q3SLX9', 'Q8DTM1', 'B8ERL9', 'Q84CG1', 'Q82WA8', 'Q824A4', 'B8J3V0', 'Q7VL09', 'B0U6J0', 'Q72LL6', 'B0UKC8', 'Q6AL81', 'B2SWS7', 'Q68VS3', 'P26505', 'Q5SHW0', 'P57202', 'Q5L6M0', 'Q56114', 'Q55128', 'Q3Z8H5', 'Q4UJV4', 'A2SD53', 'Q73KM3', 'P55683', 'Q318P3', 'Q2JS04', 'B5F073', 'Q8PQD8', 'Q7W9I4', 'Q1RIV2', 'Q2W3L2', 'Q1MR87', 'B1Y500', 'A1KVF6', 'Q08432', 'Q06191', 'P37419', 'Q7VQW9', 'A0QHJ9', 'A1VUJ6', 'Q04UL5', 'Q02636', 'P9WQ91', 'P9WQ90', 'Q89AK6', 'B2FLM5', 'B4ESU4', 'Q5NZF5', 'P97084', 'P96847', 'P95468', 'B2IH50', 'P77806', 'Q7WDY3', 'B7ULX3', 'P72173', 'A6L2V8', 'Q5F6R6', 'Q1LT68', 'P58661', 'Q2Y9Y8', 'A9HJ57', 'P44425', 'P43336', 'B4RFX5', 'P36692', 'P28735', 'P23034', 'P16524', 'P0A961', 'P0A960', 'P0A959', 'P09053', 'A9VG56', 'P00509', 'Q6D3C0', 'O85746', 'B9KDN6', 'A9BV10', 'O52815', 'O25320', 'O07587', 'A7GSE1', 'C6C2Z3', 'C0QFJ4', 'B8CX89', 'B3PI88', 'B3Q8Z5', 'Q9K625', 'B2A250', 'Q2NUJ6', 'A5FZN8', 'Q3BYN0', 'Q6LPR3', 'B0CDH5', 'Q5KY23', 'A1WVM6', 'P22806', 'A8AJ11', 'A7N6R9', 'B2VBT8', 'Q9ZLN3', 'Q3Z408', 'A1ATI6', 'Q74GT3', 'B6ESC6', 'Q81I05', 'A1K6Q1', 'Q119K2', 'P9WPZ4', 'A1AQT1', 'B5YFU5', 'A1TTV0', 'A7MX30', 'Q9PK04', 'B8F713', 'Q82DR2', 'Q06965', 'Q84I52', 'O86459', 'Q8DJ97', 'C4ZG66', 'P18079', 'Q04512', 'Q3JMZ7', 'P44422', 'A4G5N9', 'A3PMF8', 'Q1I3N7', 'Q058A6', 'P43089', 'B7VH15', 'B1XL23']

    # pf04055 kmedoids size 150
    # custom_subset = ['Q7NFJ9', 'Q51676', 'A0A384LP51', 'Q057Q7', 'Q72DS4', 'Q89AK5', 'A6LUQ4', 'Q9ZLH0', 'A6H1N2', 'Q30XT6', 'Q38HX2', 'Q2GG00', 'Q92CY2', 'Q74CF3', 'B3EGT4', 'Q30Y73', 'Q7NIT2', 'Q2JRI4', 'P07748', 'Q5SGZ3', 'B8J0C6', 'B6JNP9', 'B4SHB4', 'Q5HSB2', 'B4SA62', 'Q2RSY6', 'C1DEW5', 'P36569', 'Q2NYZ2', 'Q72DM4', 'Q7NCE3', 'Q9X758', 'A0LV48', 'B8ENI9', 'Q44634', 'Q5YZ59', 'B1HQE6', 'P0A1E2', 'B6IQ36', 'Q728P5', 'Q3AGB7', 'A7NPY6', 'Q057Q1', 'P57357', 'B2RMI0', 'A9BGV7', 'A9A0B5', 'Q6F9I9', 'A6LD84', 'A8L6D8', 'A0Q2E1', 'B8GXM4', 'A7I414', 'O31677', 'A0L509', 'B2GLQ7', 'A7GH77', 'P57373', 'Q609V2', 'A1AWE7', 'O33506', 'Q89AK8', 'B9J9I5', 'Q7U3H2', 'A9NEU7', 'A9CF16', 'Q81G67', 'B1IHH9', 'Q8DLC2', 'Q18CP3', 'Q5P7B0', 'P32131', 'Q8DII8', 'A7I1A0', 'B1X7B1', 'A8GBC5', 'Q84F14', 'A9EPV3', 'B2FUS8', 'A8G2A6', 'A8Z642', 'B3CQP5', 'B3E599', 'A9KG71', 'C1C6B0', 'Q0TTH1', 'A5VUG1', 'A9FD89', 'P72811', 'B1ZQZ5', 'Q8KD71', 'Q9K0Q5', 'B3DYX1', 'Q53U14', 'Q8K9D9', 'Q55373', 'B4ESU3', 'Q0BSW6', 'Q7VFE6', 'Q9EYN8', 'Q55914', 'Q02550', 'Q31PU7', 'Q8KFK8', 'Q8D1Y5', 'P9WJ78', 'B7VAB6', 'P45097', 'P95651', 'A3PJL9', 'A5VXF1', 'Q8EVK0', 'Q0RDQ8', 'P75794', 'P74132', 'P73667', 'Q8K9P5', 'A2SLX7', 'Q057G5', 'Q97L63', 'P55477', 'Q56211', 'P71517', 'B8D983', 'Q8XMQ3', 'Q97D55', 'Q49573', 'Q3JCN4', 'Q3AVP9', 'Q9L3B0', 'O67886', 'P59749', 'B6JNS0', 'B5ECN1', 'Q899M1', 'Q3ASD4', 'Q2RT68', 'Q9S498', 'A8Z609', 'P30140', 'B7ID25', 'Q3YRT2', 'B1XX30', 'A9KK15', 'E3PRJ8', 'Q8YR77', 'B2A3C0', 'B0SGV6', 'Q2NAD7', 'A5G3Q7']

    METHOD: str = "weighted"
    mat = np.load(_MAT)
    Z = linkage(mat, method=METHOD)

    score_lst = []
    cov_lst = []
    gsi_lst = []
    simscore_lst = []

    for i in range(1):
        sol = random_solution(len(head), 100)
        sol_subset = [list(ac_to_ec.keys())[i] for i in range(len(sol)) if sol[i] == 1]

        sol_binary = []

        for item in ac_to_ec.keys():

            if item in custom_subset:
                sol_binary.append(1)
            else:
                sol_binary.append(0)

        # sol_binary = get_n(Z, 115, mat)
        sol_subset = [list(ac_to_ec.keys())[i] for i in range(len(sol_binary)) if sol_binary[i] == 1]


        score, sim_list = score_sol(mat, sol_binary)
        # score, sim_list = score_sol(mat, sol)
        score_lst.append(score)
        simscore_lst.append(np.mean(sim_list))

        # ecs = get_ec_subset(sol_subset, ac_to_ec)
        ecs = get_ec_subset(custom_subset, ac_to_ec)

        cov = len(ecs.keys()) / float(len(ec_to_ac.keys()))

        cov_lst.append(cov)

        # gs_dict = gini_simpson_dict(sol_subset, ac_to_ec)
        gs_dict = gini_simpson_dict(custom_subset, ac_to_ec)

        gsi = gini_simpson_value(gs_dict)

        gsi_lst.append(gsi)

    print(f"Score:{np.mean(score_lst)}±{np.std(score_lst)}\n\
            Cov: {np.mean(cov_lst)}±{np.std(cov_lst)}\n\
            GSI: {np.mean(gsi_lst)}±{np.std(gsi_lst)}\n\
            SeqSim: {np.mean(simscore_lst)}±{np.std(simscore_lst)}")

    score, sim_list = score_sol(mat, sol_binary)

    gs_dict = gini_simpson_dict(list(ac_to_ec.keys()), ac_to_ec)
    gsi = gini_simpson_value(gs_dict)
    full_ec_lst = set(sum([ec_to_ac[ec] for ec in ec_to_ac.keys() if '.-' not in ec],[]))

    print(f"Num Seq:{len(ac_to_ec.keys())}\nGSI: {gsi}\nNum Ann:{len(ec_to_ac.keys())}")

    gsi_dict1 = gs_dict
    # gsi_dict2 = gini_simpson_dict(sol_subset, ac_to_ec)
    gsi_dict2 = gini_simpson_dict(custom_subset, ac_to_ec)

    gsi2 = gini_simpson_value(gsi_dict2)

    new_dict = defaultdict(float)
    diff_lst = []
    prop_lst = []
    gsi_vals1 = []
    gsi_vals2 = []

    for sig in gsi_dict1.keys():
        if sig not in gsi_dict2.keys():
            gsi_dict2[sig] = 0
        gsi_vals1.append(gsi_dict1[sig])
        gsi_vals2.append(gsi_dict2[sig])

    for sig in gsi_dict1.keys():
        diff = gsi_dict2[sig] - gsi_dict1[sig]
        diff_lst.append(diff)

        new_dict[sig] = diff
        sig_prop = (len(ec_to_ac[sig]) / len(ac_to_ec.keys()))**2
        prop_lst.append(sig_prop)
        print(f"{sig}\t{diff}\t{sig_prop}")


    # plt.hist(diff_lst, bins=np.linspace(-0.25, 0.2, num=100))

    # plt.ylim(-0.1, 0.1)
    # plt.xlim(0, 1)
    # sns.regplot(prop_lst, diff_lst, truncate=False)


    temp_prop_lst = deepcopy(prop_lst)
    print(prop_lst)
    print(gsi_vals1, gsi_vals2)
    prop_lst, gsi_vals1 = (list(t) for t in zip(*sorted(zip(prop_lst, gsi_vals1))))
    # print(prop_lst, gsi_vals1)
    prop_lst, gsi_vals2 = (list(t) for t in zip(*sorted(zip(temp_prop_lst, gsi_vals2))))
    # print(prop_lst, gsi_vals2)

    print(gsi_vals1, gsi_vals2)

    #
    # prop_lst, gsi_vals1, gsi_vals2 = map(list, zip(*sorted(zip(prop_lst, gsi_vals1, gsi_vals2), reverse=True)))
    #

    fig, ax = plt.subplots()
    plt.title(r"Change in relative abundance for $\bf{ADH}$ dataset", fontsize=16)

    prop_lst = list(range(len(gsi_vals1)))
    ax.scatter(prop_lst, gsi_vals1, label="Superset", alpha=0.8)
    ax.scatter(prop_lst, gsi_vals2, label="Subset", alpha=0.8)
    # plt.xlabel("$p_i$", fontsize=14)
    plt.ylabel("$p^2_i$", fontsize=14)
    plt.text(0.025, 0.75, "$\overline{p^2_i}$, $i \in (U\cup Z)$=" + f"{round(np.mean(gsi_vals1), 4)}", ha='left', va='top', transform=ax.transAxes, fontsize=13)
    plt.text(0.025, 0.65, "$\overline{p^2_i}$, $i \in U$=" + f"{round(np.mean(gsi_vals2), 4)}", ha='left', va='top', transform=ax.transAxes, fontsize=13)

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

    plt.legend(fontsize=12, loc='upper left')
    # print(prop_lst, gsi_vals1)
    plt.fill_between(prop_lst, gsi_vals1, gsi_vals2, step='mid', color='grey', alpha='0.3')

    print(len(gsi_vals1), len(gsi_vals2))
    plt.show()
    # plt.savefig("PF04055_deltaGSI.pdf", dpi=500)



if __name__ == '__main__':
    main()
