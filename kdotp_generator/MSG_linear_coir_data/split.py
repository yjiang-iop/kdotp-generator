import numpy as np
import os

data = np.load('MSG_linear_coir.npy',allow_pickle=True)

for msg_dict in data:
    msg_num = msg_dict[0]['msg_num']
    print(msg_num)
    np.save(msg_num+'.npy', msg_dict)
