try:
    import colorcet as cc
except:
    print "Cannot convert from colorcet w/o colorcet"
    import sys
    sys.exit()

inverse = {}
for k,v in cc.cm.items():
    if not k[-2:] == "_r":
        inverse[v] = inverse.get(v, [])
        inverse[v].insert(0,k)
all_cms = {',  '.join(reversed(v)):k for (k,v) in inverse.items()}.items()
all_cms.sort()

def dump_cmap(mpl_cmap):
    nm = mpl_cmap.name
    with open("%s.rgb" % nm,"w") as f:
        f.write("# Converted from colorcet\n")
        f.write("#\n")
        f.write("# number of colors in table\n")
        f.write("#ncolors = %i\n" % mpl_cmap.N)
        f.write("#\n")
        f.write("#  r   g   b\n")
        for i in range(mpl_cmap.N):
            a = float(i)/float(mpl_cmap.N-1)
            r,g,b,a = [int(x*255) for x in mpl_cmap(a)]
            f.write(" %3s %3s %3s\n" % (r, g,b))
    print "Wrote %s" % nm
    
for cmap in all_cms:
    dump_cmap(cmap[1])

