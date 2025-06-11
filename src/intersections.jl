function intersections(G::FiniteCoxeterGroup)
    H = hecke(G,Pol(:q)) # Hecke algebra object
    H_ct = CharTable(H).irr  # Hecke algebra character table array
    W_ct = map(x->x(1),H_ct) # Weyl group character table array

    ucl = UnipotentClasses(G) # Unipotent classes object
    gt = GreenTable(ucl;classes=true) # Green functions object
    gt_vals = Pol{Rational{Int64}}.(Pol.(gt.scalar)) # Green functions values array

    #F_mat = drinfeld_double(G).fourierMat # Fourier matrix

    uch = UnipotentCharacters(G) # Unipotent characters object
    uval = UnipotentValues(ucl;classes=true) # Table of unipotent characters on unipotent classes
    uval_ct = Pol{Rational{Int64}}.(Pol.(uval.scalar)) # Aforementioned table as an array

    # Determine which unipotent characters are principal (happens iff fake degree > 0)
    # Fake degree := substitute q=1 in degree polynomial
    fake_degs = map(p->p(1),degrees(uch))[:,1]
    principal_uch_indices = findall(d->d>0,fake_degs)
    # Remove char table rows of non-principal characters
    uval_ct = uval_ct[principal_uch_indices,begin:end]

    # Compute intersection numbers |BwB∩C|
    # (We follow Example 3.9 in [Geck 2011, Some applications of CHEVIE])
    xt = XTable(ucl;classes=true) # Contains centraliser sizes of unipotent classes
    intersections_unscaled = transpose(H_ct)*uval_ct 
    # Rescale each entry by size of Borel and centraliser sizes
    borel_size = Pol(:q)^(length(roots(G))/2)*(Pol(:q)-1)^rank(G)
    centraliser_sizes = Pol{Rational{Int64}}.(Pol.(xt.centClass))
    intersections = Pol{Rational{Int64}}.((borel_size .* intersections_unscaled) * Diagonal(map(p->1//p,centraliser_sizes)))

    # Labels for table
    rational_unipotent_classes_TeX_names = map(label -> name(TeX(rio();class=label[2]),ucl.classes[label[1]]),xt.classes)
    rational_unipotent_classes_names = fromTeX.(Ref(rio()),rational_unipotent_classes_TeX_names)
    weyl_classes_TeX_names = charnames(uch;TeX=true)[principal_uch_indices,:]
    weyl_classes_names = fromTeX.(Ref(rio()),weyl_classes_TeX_names)

    # Formatting table
    repr_data = xrepr.(Ref(rio()),CycPol.(intersections))
    println("")
    println("Rows labelled by φ∈Irr(W) <-> C∈ConjCl(W)")
    println("Columns labelled by unipotent classes C⊆G(Fq)")
    println("")
    showtable(repr_data;col_labels=rational_unipotent_classes_names,row_labels=weyl_classes_names,rows_label="|BwB∩C|")
end