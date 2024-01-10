
# moving average function
def slow_moving_average2 (signal, win):
    assert win % 2 != 0
    new_sig = []
    for i in range(len(signal)):
        neighbors = np.asarray(signal[max(0, i-win/2):min(i+win/2+1, len(signal))])
        neighbors = neighbors[~np.isnan(neighbors)]
        new_sig.append(np.mean(neighbors))
    return new_sig


#### set parameters
# set binning resolution
i = 20
bin_size = int(0.5*(10**6) / i) # binsize (unit of bp)
blur_win = int(4*i + 1) # sliding window (unit of bin)



#### binning the target data
    for name in name_ID_pos:
        ID_pos = name_ID_pos[name]
        ID_value = name_ID_value[name]
        for ID in ID_pos:
            binID = int(ID_pos[ID]) / int(bin_size)
            value = ID_value[ID]
            if np.isnan(value):
                continue
            if name not in name_binID_mean:
                name_binID_mean[name] = {}
            if binID not in name_binID_mean[name]:
                name_binID_mean[name][binID] = 0.0
            name_binID_mean[name][binID] += value
            if name not in name_binID_count:
                name_binID_count[name] = {}
            if binID not in name_binID_count[name]:
                name_binID_count[name][binID] = 0
            name_binID_count[name][binID] += 1

    for name in name_binID_mean:
        for binID in name_binID_mean[name]:
            name_binID_mean[name][binID] = float(name_binID_mean[name][binID]) / name_binID_count[name][binID]



#### smoothing the binned data by sliding window average
    name_sig = {}
    for name in names:
        binID_mean = name_binID_mean[name]
        sig = []
        for binID in range(len(gband_img)):
            #size = binID_size[binID]
            size = 1
            try:
                sig += [binID_mean[binID] for k in range(size)]
            except:
                sig += [np.nan for k in range(size)]
        if name == 'eigen':
            sig = statis.slow_moving_average2(sig, int(4*i/5.0 + 1))
        else:
            sig = statis.slow_moving_average2(sig, blur_win)
        name_sig[name] = sig
