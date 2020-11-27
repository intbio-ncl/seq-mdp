
from collections import defaultdict

_FILE = "./ipr/PF00155_trembl_ipr.tsv"


def read_file():

    fr = open(_FILE, 'r')
    res_dict = defaultdict(list)

    for line in fr:
        line = line.strip()
        temp_arr = line.split('\t')

        acc = temp_arr[0]

        if len(temp_arr) != 13:
            if temp_arr[3] == 'PANTHER':
                sig = temp_arr[4]
                res_dict[acc].append(sig)

            continue

        sig = temp_arr[11]

        res_dict[acc].append(sig)

    fr.close()

    return res_dict


def write_results(res_dict):

    fw = open('PF00155_trembl_ipr.tab', 'w')

    fw.write("Accession\tSignatures\n")

    for key in res_dict.keys():
        sigs = list(set(res_dict[key]))

        if sigs == []:
            fw.write(f"{key}\t-\n")
        else:
            sigs_str = '; '.join(sigs)
            fw.write(f"{key}\t{sigs_str}\n")


def main():

    res_dict = read_file()
    write_results(res_dict)


if __name__ == '__main__':
    main()
