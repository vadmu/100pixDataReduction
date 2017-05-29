from __future__ import division

import numpy as np
from scipy.signal import fftconvolve
from scipy.stats import scoreatpercentile, mode
from collections import deque


import peakutils
from pylab import * 
import math
from scipy import interpolate

## Finds peak positions based on parameters peak_nr (number of peaks), md (minimum distance between peaks), th(threshold of noise, 1:max, 0:min), cut (offset in the spectra):
def findPeaks(spectrum, peak_nr, md, th,cut):
	spectrum_norm=spectrum.astype(float)/max(spectrum.astype(float))*10000
	indexes= peakutils.indexes(spectrum_norm[cut:-1], thres=th, min_dist=md)
	I=np.zeros(peak_nr+1,dtype=np.int)
	I[1:peak_nr+1]=indexes[-peak_nr:]+cut
	return I
## Finds linear relation between the peak position in every pixel compared to a reference pixel:
def relatePix(I,IR):
	m,b=polyfit(I, IR, 1)
	return m,b
## Use the relation above to interpolate and align the pixels:
def align(pix, spectrum, num_channels, m, b):
	if m[pix]!= 0. or b[pix]!=0.:
		x=np.arange(0,num_channels)
		y=(m[pix]*x+b[pix]).astype(int)
		data_interp=np.interp(x, y, spectrum)
	else:
		data_interp=np.zeros((num_channels), dtype=np.int)
	return data_interp





#import matplotlib.pyplot as plt

# peak detection library:
# https://github.com/zmzhang/libPeak
#
# 1. unzipped version already works with python
# 2. to make a compiled cython version one needs to install cython
# 3. run build_pyx.py and then copy files to the folder with functions.py
# 4. in case of Ubuntu 14.04 LTS, python 2.7.6, GCC 4.8.4, numpy 1.11.3, scipy 0.18.1
#    I had to change a bit _peak_detection.pyx file and do step 3 again:
#    line 102: np.int_t --> np.int32_t 
#    line 104: np.long --> np.int  
#    line 186: comment print 

try:
    from _peak_detection import _ridge_detection, _peaks_position
    use_cython = True
except:
    use_cython = False
    
    
def mexican_hat(points, a):
    A = 2 / (np.sqrt(3 * a) * (np.pi ** 0.25))
    wsq = a ** 2
    vec = np.arange(0, points) - (points - 1.0) / 2
    tsq = vec ** 2
    mod = (1 - tsq / wsq)
    gauss = np.exp(-tsq / (2 * wsq))
    total = A * mod * gauss
    return total
 
 
def cwt(data, wavelet, widths):
    output = np.zeros([len(widths), len(data)])
    for ind, width in enumerate(widths):
        wavelet_data = wavelet(min(10 * width, len(data)), width)
        output[ind, :] = fftconvolve(data, wavelet_data,
                                     mode='same')
    return output
 
 
def local_extreme(data, comparator,
                  axis=0, order=1, mode='clip'):
    if (int(order) != order) or (order < 1):
        raise ValueError('Order must be an int >= 1')
    datalen = data.shape[axis]
    locs = np.arange(0, datalen)
    results = np.ones(data.shape, dtype=bool)
    main = data.take(locs, axis=axis, mode=mode)
    for shift in range(1, order + 1):
        plus = data.take(locs + shift, axis=axis, mode=mode)
        minus = data.take(locs - shift, axis=axis, mode=mode)
        results &= comparator(main, plus)
        results &= comparator(main, minus)
    return results
 
 
def ridge_detection(local_max, row_best, col, n_rows, n_cols, minus=True, plus=True):
    cols = deque()
    rows = deque()
    cols.append(col)
    rows.append(row_best)
    col_plus = col
    col_minus = col
    for i in range(1, n_rows):
        row_plus = row_best + i
        row_minus = row_best - i
        segment_plus = 1
        segment_minus = 1
        if minus and row_minus > 0 and segment_minus < col_minus < n_cols - segment_minus - 1:
            if local_max[row_minus, col_minus + 1]:
                col_minus += 1
            elif local_max[row_minus, col_minus - 1]:
                col_minus -= 1
            elif local_max[row_minus, col_minus]:
                col_minus = col_minus
            else:
                col_minus = -1
            if col_minus != -1:
                rows.appendleft(row_minus)
                cols.appendleft(col_minus)
        if plus and row_plus < n_rows and segment_plus < col_plus < n_cols - segment_plus - 1:
            if local_max[row_plus, col_plus + 1]:
                col_plus += 1
            elif local_max[row_plus, col_plus - 1]:
                col_plus -= 1
            elif local_max[row_plus, col_plus]:
                col_plus = col_plus
            else:
                col_plus = -1
            if col_plus != -1:
                rows.append(row_plus)
                cols.append(col_plus)
        if (minus and False == plus and col_minus == -1) or \
                (False == minus and True == plus and col_plus == -1) or \
                (True == minus and True == plus and col_plus == -1 and col_minus == -1):
            break
    return rows, cols
 
 
def peaks_position(vec, ridges, cwt2d, wnd=2):
    n_cols = cwt2d.shape[1]
    negs = cwt2d < 0
    local_minus = local_extreme(cwt2d, np.less, axis=1, order=1)

    zero_crossing = np.abs(np.diff(np.sign(cwt2d))) / 2
    # # figure(figsize=(12, 3))
    # imshow(zero_crossing, cmap=cmap_black)
    # ylabel("Scales")
    # # figure(figsize=(12, 3))
    # imshow(local_minus, cmap=cmap_blue)
    # ylabel("Scales")

    negs |= local_minus
    negs[:, [0, n_cols - 1]] = True
    ridges_select = []
    peaks = []
 
    for ridge in ridges:
        inds = np.where(cwt2d[ridge[0, :], ridge[1, :]] > 0)[0]
        if len(inds) > 0:
            col = int(mode(ridge[1, inds])[0][0])
            rows = ridge[0, :][(ridge[1, :] == col)]
            row = rows[0]
            cols_start = max(col - np.where(negs[row, 0:col][::-1])[0][0], 0)
            cols_end = min(col + np.where(negs[row, col:n_cols])[0][0], n_cols)
            # print col, row, cols_start, cols_end
            if cols_end > cols_start:
                inds = range(cols_start, cols_end)
                peaks.append(inds[np.argmax(vec[inds])])
                ridges_select.append(ridge)
        elif ridge.shape[1] > 2: # local wavelet coefficients < 0
            cols_accurate = ridge[1, 0:ridge.shape[1] / 2]
            cols_start = max(np.min(cols_accurate) - 3, 0)
            cols_end = min(np.max(cols_accurate) + 4, n_cols - 1)
            inds = range(cols_start, cols_end)
            if len(inds) > 0:
                peaks.append(inds[np.argmax(vec[inds])])
                ridges_select.append(ridge)
    # print peaks
    ridges_refine = []
    peaks_refine = []
    ridges_len = np.array([ridge.shape[1] for ridge in ridges_select])
    # print zip(peaks, ridges_len)
    for peak in np.unique(peaks):
        inds = np.where(peaks == peak)[0]
        ridge = ridges_select[inds[np.argmax(ridges_len[inds])]]
        inds = np.clip(range(peak - wnd, peak + wnd + 1), 0, len(vec) - 1)
        inds = np.delete(inds, np.where(inds == peak))
        if np.all(vec[peak] > vec[inds]):
            ridges_refine.append(ridge)
            peaks_refine.append(peak)
    return peaks_refine, ridges_refine
 
 
def ridges_detection(cwt2d, vec):
    n_rows = cwt2d.shape[0]
    n_cols = cwt2d.shape[1]
    local_max = local_extreme(cwt2d, np.greater, axis=1, order=1)
    ridges = []
    rows_init = np.array(range(1, 6))
    cols_small_peaks = np.where(np.sum(local_max[rows_init, :], axis=0) > 0)[0]
    for col in cols_small_peaks:
        best_rows = rows_init[np.where(local_max[rows_init, col])[0]]
        if use_cython:
            rows, cols = _ridge_detection(local_max, best_rows[0], col, n_rows, n_cols, True, True)
        else:
            rows, cols = ridge_detection(local_max, best_rows[0], col, n_rows, n_cols, True, True)        
        staightness = 1 - float(sum(abs(np.diff(cols)))) / float(len(cols))
        if len(rows) >= 2 and \
            staightness > 0.2 and  \
            not(
            len(ridges) > 0 and
            rows[0] == ridges[-1][0, 0] and
            rows[-1] == ridges[-1][0, -1] and
            cols[0] == ridges[-1][1, 0] and
            cols[-1] == ridges[-1][1, -1] and
            len(rows) == ridges[-1].shape[1]
        ):
            ridges.append(np.array([rows, cols], dtype=np.int32))

    # figure(figsize=(12, 3))
    # imshow(cwt2d)
    # ylabel("Scales")
#    figure()
#    plot(cwt2d[3,:])
#     figure(figsize=(9, 4))
#     imshow(local_max, cmap=cmap_red)
#     ylabel("Scales")

    return ridges
 
 
def signal_noise_ratio(cwt2d, ridges, peaks):
    n_cols = cwt2d.shape[1]
    row_one = cwt2d[0, :]
    row_one_del = np.delete(row_one, np.where(abs(row_one) < 10e-5))
    t = 3 * np.median(np.abs(row_one_del - np.median(row_one_del))) / 0.67
    row_one[row_one > t] = t
    row_one[row_one < -t] = -t
    noises = np.zeros(len(peaks))
    signals = np.zeros(len(peaks))
    for ind, val in enumerate(peaks):
        hf_window = ridges[ind].shape[1] * 1
        window = range(int(max([val - hf_window, 0])), int(min([val + hf_window, n_cols])))
        noises[ind] = scoreatpercentile(np.abs(row_one[window]), per=90)
        signals[ind] = np.max(cwt2d[ridges[ind][0, :], ridges[ind][1, :]])
    sig = [1 if s > 0 and n >= 0 else - 1 for s, n in zip(signals, noises)]
 
    # print zip(peaks, signals, noises)
    # figure()
    # plot(row_one, label='scale = 1')
    # # plot(cwt2d[1, :], label='scale = 2')
    # plot(cwt2d[3, :], label='scale = 4')
    # # plot(cwt2d[5, :], label='scale = 6')
    # legend()
    return np.sqrt(np.abs((signals + np.finfo(float).eps) / (noises + np.finfo(float).eps))) * sig, signals
 
 
def peaks_detection(vec, scales, min_snr=3):
    cwt2d = cwt(vec, mexican_hat, scales)
    ridges = ridges_detection(cwt2d, vec)
    # print ridges
    if use_cython:
        peaks, ridges = _peaks_position(vec, ridges, cwt2d)
    else:
        peaks, ridges = peaks_position(vec, ridges, cwt2d)
    
    # print ridges
    # print peaks
    snr, signals = signal_noise_ratio(cwt2d, ridges, peaks)
    # print zip(peaks, snr)
    # peaks_refine = [peak for i, peak in enumerate(peaks) if snr[i] >= min_snr]
    # signals_refine = [signal for i, signal in enumerate(signals) if snr[i] >= min_snr]
    # print peaks_refine
    peaks_refine = [peak for i, peak in enumerate(peaks) if signals[i] >= min_snr]
    signals_refine = [signal for i, signal in enumerate(signals) if signals[i] >= min_snr]
    return peaks_refine, signals_refine


#x = np.arange(0, 10, 0.05)
#y = np.sin(2.*3.14*x)
#peak_ind, sig = peaks_detection(y, np.arange(1, 30), 0.1)
#print y[peak_ind], sig
#plt.plot(x, y, "r-")
#plt.plot(x[peak_ind], y[peak_ind], "ko")
#plt.show()
#  

