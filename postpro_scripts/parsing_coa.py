
def parsing_coa (fname):
    plate_info_list = {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        cols = line.strip().split(',')
        plate_name, amount, well = cols[5], cols[-6], cols[-1]
        #cols = line.strip().split('\t')
        #plate_name, amount, well = cols[0], cols[-6], cols[5]
        plate_name = "_".join(plate_name.strip('"').split("_")[:-1])
        #amount = float(amount)
        amount = float(amount.strip('"')[:-6])
        well = well.strip('"')
        if plate_name not in plate_info_list:
            plate_info_list[plate_name] = []
        plate_info_list[plate_name].append((well, amount, amount*10))

    return plate_info_list

plate_info_list = parsing_coa("coa.csv")

for plate_name in plate_info_list:
    f = open(plate_name + '_info.txt', 'w')
    print >> f, plate_name
    print >> f, "Well\tAmount\tAddFor100uM"
    info_list = plate_info_list[plate_name]
    for info in info_list:
        print >> f, "%s\t%s nmole\t%s ul" % info
    f.close()
