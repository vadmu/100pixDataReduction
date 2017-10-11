import numpy as np

def puindexes(y, thres=0.3, min_dist=1):
    """Peak detection routine.
    Finds the numeric index of the peaks in *y* by taking its first order difference. By using
    *thres* and *min_dist* parameters, it is possible to reduce the number of
    detected peaks. *y* must be signed.
    Parameters
    ----------
    y : ndarray (signed)
        1D amplitude data to search for peaks.
    thres : float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    min_dist : int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.
    Returns
    -------
    ndarray
        Array containing the numeric indexes of the peaks that were detected
    """
    if isinstance(y, np.ndarray) and np.issubdtype(y.dtype, np.unsignedinteger):
        raise ValueError("y must be signed")

    thres = thres * (np.max(y) - np.min(y)) + np.min(y)
    min_dist = int(min_dist)

    # compute first order difference
    dy = np.diff(y)

    # propagate left and right values successively to fill all plateau pixels (0-value)
    zeros,=np.where(dy == 0)
    
    while len(zeros):
        # add pixels 2 by 2 to propagate left and right value onto the zero-value pixel
        zerosr = np.hstack([dy[1:], 0.])
        zerosl = np.hstack([0., dy[:-1]])

        # replace 0 with right value if non zero
        dy[zeros]=zerosr[zeros]
        zeros,=np.where(dy == 0)

        # replace 0 with left value if non zero
        dy[zeros]=zerosl[zeros]
        zeros,=np.where(dy == 0)

    # find the peaks by using the first order difference
    peaks = np.where((np.hstack([dy, 0.]) < 0.)
                     & (np.hstack([0., dy]) > 0.)
                     & (y > thres))[0]

    if peaks.size > 1 and min_dist > 1:
        highest = peaks[np.argsort(y[peaks])][::-1]
        rem = np.ones(y.size, dtype=bool)
        rem[peaks] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peak] = False

        peaks = np.arange(y.size)[~rem]
    return peaks

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth    
    
# Finds peak positions based on parameters peak_nr (number of peaks), md (minimum distance between peaks), th(threshold of noise, 1:max, 0:min), cut (offset in the spectra):
def findPeaks(spectrum, peak_nr, md, th, cut):
#    spectrum_norm=spectrum.astype(float)/max(spectrum.astype(float))*10000
    indexes= puindexes(smooth(spectrum[cut:], 5), thres=th, min_dist=md)
#	I=np.zeros(peak_nr+1,dtype=np.int)
#	I[1:peak_nr+1]=indexes[-peak_nr:]+cut
    highest = np.argsort(spectrum[indexes + cut])[:peak_nr]
    I = indexes[sorted(highest)] + cut
    return I


    
# Finds linear relation between the peak position in every pixel compared to a reference pixel:
def relatePix(I,IR):
    m, b=np.polyfit(I, IR, 1)
    return m, b
    
# Use the relation above to interpolate and align the pixels:
def align(pix, spectrum, num_channels, m, b):
    if m[pix]!= 0. or b[pix]!=0.:
        x=np.arange(0,num_channels)
        y=m[pix]*x+b[pix]
        data_interp=np.interp(x, y, spectrum)
    else:
        data_interp=np.zeros((num_channels), dtype=np.int)
    return data_interp
    
 