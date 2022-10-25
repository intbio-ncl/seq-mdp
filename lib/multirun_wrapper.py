
import subprocess
import numpy as np

_I = 50

_OUT = "./PF00155_ts_runs_csn.tab"


def main():

    cov_lst = []
    gsi_lst = []
    fw = open(_OUT, 'w')
    fw.write("EC Coverage\tEC GSI\n")

    for i in range(_I):
        print(f"RUN NUMBER {i+1}")
        fw = open(_OUT, 'a')
        subprocess.run("python3 main.py -hd PF00155_csn_headings.json -d PF00155_csn_mat.npy -k 100  -s ts -a PF00155_annotation.tab", shell=True)
        fr = open("ts_vals.txt", "r")
        out = list(fr)
        cov = float(out[0].split(' ')[-1].strip())
        gsi = float(out[1].split(' ')[-1].strip())
        cov_lst.append(cov)
        gsi_lst.append(gsi)
        fw.write(f"{cov}\t{gsi}\n")
        fw.close()
        fr.close()
    print(cov_lst)
    print(gsi_lst)
    print(f"EC Coverage: {np.mean(cov_lst)}±{np.std(cov_lst)}")
    print(f"EC GSI: {np.mean(gsi_lst)}±{np.std(gsi_lst)}")

    fw.close()


if __name__ == '__main__':
    main()
