1. copy the hypersurface files from CCAKE to ./BSQ-FO-ideal/input/hs_.../ev${i}
2. modify your dfinput.dat 
3. run the following from the top level directory: 
 bash submit_batches ./out/hs_myfreezeout_folder/ hs_myfreezeout_folder N
N here is the number of events you are attempting to freezeout.

Note: the hs_myfreezeout_folder requires the following format:
1. ev0.......evN for the subfolders which includes the hypersurfaces
2. rename your hypersurface_i.dat to sbvfreezeout0.dat