import requests
import pymol

def getImage(pdb_name, dir):
    # a generic method to get content from a url
    url = f"https://cdn.rcsb.org/images/structures/{pdb_name}_chain-A.jpeg"
    response = requests.get(url)
    if response.status_code == 200:  # if success
        with open(f"{dir}/{pdb_name}.jpeg", 'wb') as f:  # create a file
            f.write(response.content)

def makeImage(pdb_list):
    pymol.finish_launching(['pymol', '-cq']) # lunach a pymol thread in quiet mode
    for pdf_file, pdb_id in pdb_list:
        pymol.cmd.load(pdf_file, pdb_id)
        pymol.cmd.disable("all")
        pymol.cmd.enable(pdb_id)
        #print(pymol.cmd.get_names())
        pymol.cmd.hide('all')
        pymol.cmd.show('cartoon')
        pymol.cmd.set('ray_opaque_background', 0)
        pymol.cmd.color('red', 'ss h')
        pymol.cmd.color('yellow', 'ss s')
        pymol.cmd.png(f"{pdf_file}.png")
    pymol.cmd.quit()

if __name__ == '__main__':
    makeImage([("C:/Temp/pdb/test/281_36031_1C0A_B.pdb", "281_36031_1C0A_B"),
               ("C:/Temp/pdb/test/281_36031_1C04_E.pdb", "281_36031_1C04_E")])