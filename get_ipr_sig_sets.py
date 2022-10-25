
_DATASET = "ADH"

tabu_set = ["IPR012134", "PTHR43720:SF1", "PTHR42804", "PTHR43353", "IPR016163", "PTHR43353:SF6", "PTHR11699", "PTHR11699:SF68", "IPR016160", "IPR017628", "PTHR43217:SF2", "IPR016161", "IPR023510", "IPR015590", "PTHR43111:SF1", "PTHR42991", "PTHR42991:SF2", "PTHR11699:SF228", "PTHR11699:SF222", "PTHR43217", "PTHR43570:SF4", "PTHR43570:SF16", "PTHR43720", "PTHR42986", "IPR011985", "PTHR11063", "IPR020593", "PTHR43353:SF5", "IPR010102", "PTHR11699:SF211", "IPR029510", "PTHR43866:SF4", "IPR010061", "PTHR11699:SF266", "IPR016162", "PTHR11063:SF8", "PTHR43111", "PTHR11699:SF115", "PTHR43217:SF1", "PTHR43570", "IPR000965", "PTHR43866:SF8", "IPR015657", "IPR017749", "IPR017656", "IPR017649", "PTHR42986:SF1", "PTHR43353:SF4", "IPR011264", "IPR012394", "PTHR42991:SF1"]

greedy_set = ['IPR000965', 'IPR010061', 'IPR010102', 'IPR011264', 'IPR011985', 'IPR012134', 'IPR012394', 'IPR015590', 'IPR015657', 'IPR016160', 'IPR016161', 'IPR016162', 'IPR016163', 'IPR017628', 'IPR017649', 'IPR017656', 'IPR017749', 'IPR020593', 'IPR023510', 'IPR029510', 'PTHR11063', 'PTHR11063:SF8', 'PTHR11699', 'PTHR11699:SF211', 'PTHR11699:SF222', 'PTHR11699:SF228', 'PTHR11699:SF266', 'PTHR11699:SF68', 'PTHR42804', 'PTHR42986', 'PTHR42986:SF1', 'PTHR42991', 'PTHR42991:SF1', 'PTHR42991:SF2', 'PTHR43111', 'PTHR43111:SF1', 'PTHR43217', 'PTHR43217:SF1', 'PTHR43217:SF2', 'PTHR43353', 'PTHR43353:SF4', 'PTHR43353:SF5', 'PTHR43353:SF6', 'PTHR43521', 'PTHR43521:SF2', 'PTHR43570', 'PTHR43570:SF16', 'PTHR43570:SF4', 'PTHR43720', 'PTHR43720:SF1', 'PTHR43860', 'PTHR43860:SF2', 'PTHR43866:SF4', 'PTHR43866:SF8']


intersection = list(set(tabu_set) & set(greedy_set))
ts_unique = list(set(tabu_set) - set(greedy_set))
greedy_unique = list(set(greedy_set) - set(tabu_set))

fw_int = open(f'./ipr_set_lists/{_DATASET}_intersections.txt', 'w')
fw_unique_greedy = open(f'./ipr_set_lists/{_DATASET}_unique_greedy.txt', 'w')
fw_unique_ts = open(f'./ipr_set_lists/{_DATASET}_unique_ts.txt', 'w')

[fw_int.write(f'{ipr}\n') for ipr in intersection]
[fw_unique_greedy.write(f'{ipr}\n') for ipr in greedy_unique]
[fw_unique_ts.write(f'{ipr}\n') for ipr in ts_unique]

fw_int.close()
fw_unique_ts.close()
fw_unique_greedy.close()
