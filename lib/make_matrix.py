
import numpy as np
from collections import defaultdict
import json

_INPUT = "./PF00155_trembl_mat.txt"


def read_identities(file):

    fr = open(file, 'r')
    mat = np.zeros((10_000, 10_000))
    mat_dict = defaultdict(dict)
    index_dict = defaultdict(int)

    counter = 0
    lines = list(fr)

    for line in lines:
        line = line.strip()
        temp_arr = line.split(' ')

        mat_dict[temp_arr[0]][temp_arr[1]] = temp_arr[2]
        mat_dict[temp_arr[1]][temp_arr[0]] = temp_arr[2]

    counter = 0

    headers = sorted(mat_dict)

    for key in sorted(mat_dict):
        index_dict[key] = counter
        counter += 1

    for key_i in mat_dict.keys():
        for key_j in mat_dict[key_i]:
            i = index_dict[key_i]
            j = index_dict[key_j]

            mat[i, j] = mat_dict[key_i][key_j]
            mat[j, i] = mat_dict[key_i][key_j]

    return mat, headers


def save_matrix(mat):

    np.save('PF00155_trembl_mat.npy', mat)


def save_header(headers):

    with open(f'PF00155_trembl_headings.json', 'w') as outfile:
        json.dump(list(headers), outfile)


def main():

    mat, headers = read_identities(_INPUT)

    save_matrix(mat)
    save_header(headers)


if __name__ == '__main__':
    main()
