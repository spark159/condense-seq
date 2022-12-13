import re
import load_file

#cell_choice = ['mESCs (mouse embryonic stem cells)']

cell_choice = None
specie_choice = ['Mouse (Mus musculus)']

GEO_set = set([])

#fname = 'gtrd_prc1.txt'
#fname = 'gtrd_prc2.txt'
fname = 'gtrd_foxp3.txt'
ID_field_value = load_file.read_tabular_file (fname)

for ID in ID_field_value:
    cell = ID_field_value[ID]['Cell']
    specie = ID_field_value[ID]['Specie']

    if cell_choice!=None and cell not in cell_choice:
        continue
    if specie_choice!=None and specie not in specie_choice:
        continue

    GEO_list = ID_field_value[ID]['External References'].split(',')
    GEO_list = [code.strip() for code in GEO_list]
    
    prev_head = ""
    for code in GEO_list:
        head, number = re.split('(\d+)', code)[:-1]
        if len(head) > 0:
            GEO_set.add((head, int(number)))
            if head != prev_head:
                prev_head = head
        else:
            GEO_set.add((prev_head, int(number)))


f = open(fname.rsplit('.',1)[0] + '_output.txt', 'w')
for head, number in sorted(list(GEO_set)):
    print >> f, head + str(number)
f.close()

