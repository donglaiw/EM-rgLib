from rg_util import load_aff_plane

def join_segs(block_files, join_files, thd, max_count):
    '''Create and process the region graph
    :param block_files: the region graphs from the segmentations
    :param join_files: the region graphs from the block joins
    :param thd: the minimum threshold allowed when doing a join
    :param max_count: the maxiumum number of voxels in an object
    :returns: a dictionary mapping x, y and z to a map of segmentation ID to
    global segmentation ID.
    '''
    off = 0
    doff = {}
    dcounts = {}
    rg_affs = []
    id1s = []
    id2s = []
    print 'load block..'
    for block_file in block_files:
        b = json.load(open(block_file))
        x, y, z = [b[_] for _ in "xyz"]
        counts, rg_aff, id1, id2 = \
            [np.array(b[_]) for _ in "counts", "rg_affs", "id1", "id2"]
        doff[x, y, z] = off
        w = np.where(rg_aff >= thd)
        if len(w[0] > 0):
            rg_aff, id1, id2 = rg_aff[w], id1[w], id2[w]
            id1 = id1 + off
            id2 = id2 + off
            rg_affs.append(rg_aff)
            id1s.append(id1.astype(np.uint32))
            id2s.append(id2.astype(np.uint32))
        dcounts[x, y, z] = counts
        off += len(counts)
    counts = np.zeros(off, np.uint32)
    for (x, y, z), off in doff.items():
        c = dcounts[x, y, z]
        counts[off:off+len(c)] = c
    print 'load json..'
    for join_file in join_files:
        jd = json.load(open(join_file))
        k1 = (jd["x1"], jd["y1"], jd["z1"])
        k2 = (jd["x2"], jd["y2"], jd["z2"])
        id1 = np.array(jd["id1"]) + doff[k1]
        id2 = np.array(jd["id2"]) + doff[k2]
        rg_aff = np.array(jd["rg_affs"])
        w = np.where(rg_aff >= thd)
        id1s.append(id1[w])
        id2s.append(id2[w])
        rg_affs.append(rg_aff[w])
    rg_affs = np.hstack(rg_affs)
    id1s = np.hstack(id1s)
    id2s = np.hstack(id2s)
    #
    print "Sort by descending affinity"
    #
    order = np.lexsort((id2s, id1s, 1-rg_affs))
    rg_affs, id1s, id2s = rg_affs[order], id1s[order], id2s[order]
    #
    print "Find the global mapping"

    mapping = zwatershed.zw_do_mapping(
        id1s.astype(np.uint64),
        id2s.astype(np.uint64),
        counts.astype(np.uint64),
        max_count)
    result = {}
    for x, y, z in doff:
        n_ids = len(dcounts[x, y, z])
        off = doff[x, y, z]
        m = mapping[off:off+n_ids]
        m[0] = 0
        result[x, y, z] = m
    return result


def compute_rg(plane, sega, segb, aff_plane):
    '''Join two segmentations adjacent in X,Y,Z
    
    :param sega: the segmentation on the left
    :param segb: the segmentation on the right
    :param aff_plane: the affinity value between two volumes
    :returns: a 3-tuple of highest affinity connecting segments in sega and
        segb, the ID of the segment in sega and the ID of the segment in segb.
    '''
    if plane == 'x':
        sega_plane = sega[:, :, -1]
        segb_plane = segb[:, :, 0]
    elif plane == 'y':
        sega_plane = sega[:, -1]
        segb_plane = segb[:, 0]
    elif plane == 'z':
        sega_plane = sega[-1]
        segb_plane = segb[0]
    return compute_rg_common(sega, segb, sega_plane, segb_plane, aff_plane)

def compute_rg_common(sega, segb, sega_plane, segb_plane, aff_plane):
    '''Utility function for all 3 join directions'''
    n_a = np.max(sega)
    n_b = np.max(segb)
    mask = (sega_plane !=0) & (segb_plane != 0)
    seg_plane = sega_plane.astype(np.uint64) + n_a * segb_plane.astype(np.uint64)
    seg_plane[~ mask] = 0
    all_ids = np.unique(seg_plane[mask])
    max_aff = scipy.ndimage.maximum(aff_plane, seg_plane, all_ids)
    id1 = all_ids % n_a
    id1[id1 == 0] = n_a
    id2 = all_ids / n_a
    order = np.lexsort((id2, id1, 1-max_aff))
    id1, id2, max_aff = [_[order] for _ in id1, id2, max_aff]
    return max_aff, id1, id2
