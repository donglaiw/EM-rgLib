def get_blockId(x0, x1, sz, offset=0):
    # input: start and end coordiante
    # output: block id
    # take floor
    xi0 = (x0-offset) // sz 
    # take ceil
    xi1 = (x1-offset+sz-1) // sz
    return xi0, xi1

def load_aff_plane(fn, plane, sz, x,y,z, aff_opt=0)
    if plane == 'x':
        aff_plane = load_aff(fn, x+sz[2]-1, x+sz[2], 
                             y, y+sz[1], 
                             z, z+sz[0])[:, :, 0]
    elif plane == 'y':
        aff_plane = load_aff(fn, x, x+sz[2], 
                             y + sz[1]-1, y+sz[1], 
                             z, z+sz[0])[:, 0]
    elif plane == 'z':
        aff_plane = load_aff(fn, plane, x, x+sz[2], 
                             y, y+sz[1], 
                             z+sz[0]-1, z+sz[0])[0]
    return aff_plane


def load_aff(fn, x0, x1, y0, y1, z0, z1, aff_opt):
    '''Load a volume of x, y or z affinity
    :param xyz: either "x", "y" or "z" to pick the affinity channel to load
    :param x0: the x start coordinate.
    :param x1: the x end coordinate.
    :param y0: the y start coordinate.
    :param y1: the y end coordinate.
    :param z0: the z start coordinate.
    :param z1: the z end coordinate.
    '''
    xi0, xi1 = get_blockId(x0, x1, bSize[2], bOffset[2])
    yi0, yi1 = get_blockId(y0, y1, bSize[1], bOffset[1])
    zi0, zi1 = get_blockId(z0, z1, bSize[0], bOffset[0])

    result = np.zeros((z1 - z0, y1 - y0, x1 - x0), np.uint8)
    for xi in range(xi0, xi1):
        for yi in range(yi0, yi1):
            for zi in range(zi0, zi1):
                x0b = xi * bSize[2] + bOffset[2]
                y0b = yi * bSize[1] + bOffset[1]
                z0b = zi * bSize[0] + bOffset[0]
                block = tifffile.imread(d[xyz, x0b, y0b, z0b])
                x0v = max(x0b, x0)
                y0v = max(y0b, y0)
                z0v = max(z0b, z0)
                x1v = min(x0b+bSize[2], x1)
                y1v = min(y0b+bSize[1], y1)
                z1v = min(z0b+bSize[0], z1)
                result[z0v - z0: z1v - z0,
                       y0v - y0: y1v - y0,
                       x0v - x0: x1v - x0] = \
                    block[z0v - z0b: z1v - z0b,
                          y0v - y0b: y1v - y0b,
                          x0v - x0b: x1v - x0b]
    return result.astype(np.float32) / 255
