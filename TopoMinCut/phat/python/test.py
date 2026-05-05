import phat

boundary_matrix = phat.boundary_matrix()


boundary_matrix.columns = [ (0, []),
                            (0, []),
                            (1, [0,1]),
                            (0, []),
                            (1, [1,3]),
                            (1, [0,3]),
                            (2, [2,4,5])]

print("\nThe boundary matrix has %d columns:" % len(boundary_matrix.columns))
for col in boundary_matrix.columns:
    s = "Column %d represents a cell of dimension %d." % (col.index, col.dimension)
    if (col.boundary):
        s = s + " Its boundary consists of the cells " + " ".join([str(c) for c in col.boundary])
    print(s)
print("Overall, the boundary matrix has %d entries." % len(boundary_matrix))

pairs = boundary_matrix.compute_persistence_pairs(reduction=phat.reductions.standard_reduction)

pairs.sort()

print("\nThe boundary matrix has %d columns:" % len(boundary_matrix.columns))
for col in boundary_matrix.columns:
    s = "Column %d represents a cell of dimension %d." % (col.index, col.dimension)
    if (col.boundary):
        s = s + " Its boundary consists of the cells " + " ".join([str(c) for c in col.boundary])
    print(s)
print("Overall, the boundary matrix has %d entries." % len(boundary_matrix))


print("\nThere are %d persistence pairs: " % len(pairs))
for pair in pairs:
    print("Birth: %d, Death: %d" % pair)