python ./rd_data_3D.py --SFs $(cat scalefactors.txt) --SEQ $(awk '{print $1}' sequences.txt) --g $(awk '{print $2}' sequences.txt) --d $(awk '{print $3}' sequences.txt) --pot "$1"
