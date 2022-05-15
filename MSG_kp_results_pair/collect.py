import numpy as np 
import os

Library = '../data'
msg_data = np.load(Library+"/msg_data.npy",allow_pickle=True)



def mv():
    home = os.getcwd()
    for msg in range(1651):
        msg_dict = msg_data[msg]
        msg_num = msg_dict['msg_num']
        sg = int(msg_dict['msg_num'].split('.')[0])
        msg_type = msg_dict['msg_type']
        if msg_type != 4:
            continue
        
        print(msg_num)
       #os.system('mkdir '+msg_num)

        names = ['k']#,'E','B','kE','kB', 'EB', 'epsilon']
        for name in names:
            os.system('cp ../batch_pairs/batch_%s/%s/result_%s ./%s/'%(name, msg_num, name, msg_num)) 

mv()

