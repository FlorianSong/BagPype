import urllib.request
import os 


url = "ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif"
file_name = "components.cif"
urllib.request.urlretrieve(url, file_name)

current_file_connection = open("temp", "w")
try:
    os.stat("mmcif")
except:
    os.mkdir("mmcif")
with open(file_name) as f:
    for line in f:
        if line.startswith("data"):
            current_file_connection.close()
            current_code = line[-4:].strip()
            print(current_code)
            current_file = "mmcif/" + current_code + ".cif"
            current_file_connection = open(current_file, "w")
            current_file_connection.write(line)
        else:
            current_file_connection.write(line)

os.remove("temp")
os.remove(file_name)