###############################################################################
#####                           hist2otsu                                 #####
#####                        Evan Lahr, 2021                              #####
#####         Compute an Ostu threshold using only histogram data         #####
#####         Similar to MATLAB function 'otsuthresh' (link below)        #####
#####       https://www.mathworks.com/help/images/ref/otsuthresh.html     #####
#####                                                                     #####
#####   INPUT DATA must be a 2-column array of histogram information      #####
#####       COLUMN 1 must be the LOWER BIN EDGES of the histogram         #####
#####       COLUMN 2 must be the NORMALIZED COUNTS of the histogram       #####
##### OUTPUT is the otsu threshold value (thresh min intraclass variance) #####
###############################################################################


def hist2otsu(hist):
    bins    = np.array([item[0] for item in hist])
    counts  = np.array([item[1] for item in hist])
    binstep = bins[1] - bins[0]
    midpts  = np.array([x+(binstep/2) for x in bins])
    total_weight   = sum(counts)
    least_variance = -1
    least_variance_threshold = -1
    #create an array of all possible threshold values to loop through
    thresholds = np.arange(np.min(bins), np.max(bins), binstep)
    #loop through thresholds to find the minimum intraclass class variance
    for i in thresholds:
        bg_midpts = midpts[bins <= i]
        bg_counts = counts[bins <= i]
        if sum(bg_counts) == 0:
          continue
        weight_bg = sum(bg_counts) / total_weight
        mean_bg = np.average(bg_midpts, weights=bg_counts)
        variance_bg = np.average((bg_midpts - mean_bg)**2, weights=bg_counts)
        fg_midpts = midpts[bins > i]
        fg_counts = counts[bins > i]
        if sum(fg_counts) == 0:
           continue
        weight_fg = sum(fg_counts) / total_weight
        mean_fg = np.average(fg_midpts, weights=fg_counts)
        variance_fg = np.average((fg_midpts - mean_fg)**2, weights=fg_counts)
        within_class_variance = weight_fg*variance_fg + weight_bg*variance_bg
        if least_variance == -1 or least_variance > within_class_variance:
            least_variance = within_class_variance
            least_variance_threshold = i
            #print("var:", within_class_variance, "  thresh:", i)
    return least_variance_threshold