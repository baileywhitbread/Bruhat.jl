function intersections_rational(G::FiniteCoxeterGroup)
    # Gather required data
    H = hecke(G,Pol(:q)) # Hecke algebra object
    H_ct = CharTable(H).irr  # Hecke algebra character table array
    W_ct = map(x->x(1),H_ct) # Weyl group character table array
    ucl = UnipotentClasses(G) # Unipotent classes object
    gt = GreenTable(ucl;classes=true) # Green functions object
    gt_vals = Pol{Rational{Int64}}.(Pol.(gt.scalar)) # Green functions values array
    uch = UnipotentCharacters(G) # Unipotent characters object
    uval = UnipotentValues(ucl;classes=true) # Table of unipotent characters on unipotent classes
    uval_ct = Pol{Rational{Int64}}.(Pol.(uval.scalar)) # Aforementioned table as an array
    borel_size = Pol(:q)^(length(roots(G))/2)*(Pol(:q)-1)^rank(G) # |B(Fq)| = |U(Fq)| |T(Fq)| = q^|{positive roots}| (q-1)^rank(G)

    # Determine which unipotent characters are principal (happens iff fake degree > 0)
    # Fake degree := substitute q=1 in degree polynomial
    fake_degs = map(p->p(1),degrees(uch))[:,1]
    principal_uch_indices = findall(d->d>0,fake_degs)
    uval_ct = uval_ct[principal_uch_indices,begin:end] # Removing char table rows of non-principal characters

    # Compute intersection numbers |BwB∩C| following Example 3.9 in [Geck 2011, Some applications of CHEVIE]
    intersection_numbers_unscaled = transpose(H_ct)*uval_ct 

    # Rescale each entry by size of Borel and centraliser sizes
    xt = XTable(ucl;classes=true) # Contains centraliser sizes of unipotent classes
    rational_classes_centraliser_sizes = xt.centClass
    centraliser_sizes = Pol{Rational{Int64}}.(Pol.(rational_classes_centraliser_sizes))
    intersection_numbers_scaled = Pol{Rational{Int64}}.((borel_size .* intersection_numbers_unscaled) * Diagonal(map(p->1//p,centraliser_sizes)))

    # Labels for table
    rational_unipotent_classes_TeX_names = map(label -> name(TeX(rio();class=label[2]),ucl.classes[label[1]]),xt.classes)
    rational_unipotent_classes_names = fromTeX.(Ref(rio()),rational_unipotent_classes_TeX_names)
    weyl_classes_TeX_names = classnames(G;TeX=true)
    weyl_classes_names = fromTeX.(Ref(rio()),weyl_classes_TeX_names)

    # Formatting table
    repr_data = xrepr.(Ref(rio()),CycPol.(intersection_numbers_scaled))
    println("")
    println("Rows labelled by conj classes C∈ConjCl(W)")
    println("Columns labelled by unipotent classes C⊆G(Fq)")
    println("")
    showtable(repr_data;col_labels=rational_unipotent_classes_names,row_labels=weyl_classes_names,rows_label="|BwB∩C|")
end

function intersections_geometric(G::FiniteCoxeterGroup)
    # Gather required data
    H = hecke(G,Pol(:q)) # Hecke algebra object
    H_ct = CharTable(H).irr  # Hecke algebra character table array
    W_ct = map(x->x(1),H_ct) # Weyl group character table array
    ucl = UnipotentClasses(G) # Unipotent classes object
    gt = GreenTable(ucl;classes=true) # Green functions object
    gt_vals = Pol{Rational{Int64}}.(Pol.(gt.scalar)) # Green functions values array
    uch = UnipotentCharacters(G) # Unipotent characters object
    uval = UnipotentValues(ucl;classes=true) # Table of unipotent characters on unipotent classes
    uval_ct = Pol{Rational{Int64}}.(Pol.(uval.scalar)) # Aforementioned table as an array
    borel_size = Pol(:q)^(length(roots(G))/2)*(Pol(:q)-1)^rank(G) # |B(Fq)| = |U(Fq)| |T(Fq)| = q^|{positive roots}| (q-1)^rank(G)
    xt = XTable(ucl;classes=true) # Contains information about rational unipotent classes and geometric unipotent classes

    # Determine which unipotent characters are principal (happens iff fake degree > 0)
    # Fake degree := substitute q=1 in degree polynomial
    fake_degs = map(p->p(1),degrees(uch))[:,1]
    principal_uch_indices = findall(d->d>0,fake_degs)
    uval_ct = uval_ct[principal_uch_indices,begin:end] # Removing char table rows of non-principal characters

    # Compute intersection numbers |BwB∩C| following Example 3.9 in [Geck 2011, Some applications of CHEVIE]
    intersection_numbers_unscaled = transpose(H_ct)*uval_ct 

    # Rescale each entry by size of Borel and centraliser sizes
    rational_classes_centraliser_sizes = xt.centClass
    centraliser_sizes = Pol{Rational{Int64}}.(Pol.(rational_classes_centraliser_sizes))
    intersection_numbers_scaled_unsummed = Pol{Rational{Int64}}.((borel_size .* intersection_numbers_unscaled) * Diagonal(map(p->1//p,centraliser_sizes)))

    # Sum columns corresponding to the same geometric unipotent class
    # This is where this function diverges from the intersections_rational function
    # Determine which rational classes correspond to the same geometric class
    rational_geometric_indices = xt.classes # A list of pairs [n,m] with n counting geometric orbits and m counting rational orbits inside the geometric one  
    class_ids = map(x->x[1],rational_geometric_indices) # Cuts off m
    duplicated = filter(u -> count(==(u), class_ids) > 1, unique(class_ids)) # Find all labels that occur more than once
    groups = reverse([ findall(==(u), class_ids) for u in duplicated ]) # Groups of columns which need to be summed
    
    for group in groups
        summed_col = sum(intersection_numbers_scaled_unsummed[:, group], dims=2)
        cols_before_summed_col = intersection_numbers_scaled_unsummed[:, 1:group[1]-1]
        cols_after_summed_col = intersection_numbers_scaled_unsummed[:,group[end]+1:end]
        intersection_numbers_scaled_unsummed = hcat(cols_before_summed_col, summed_col, cols_after_summed_col)
    end

    intersection_numbers_scaled_summed = intersection_numbers_scaled_unsummed # Readability

    # Labels for table
    rational_unipotent_classes_TeX_names = map(label -> name(TeX(rio();class=label[2]),ucl.classes[label[1]]),xt.classes)
    rational_unipotent_classes_names = fromTeX.(Ref(rio()),rational_unipotent_classes_TeX_names)
    weyl_classes_TeX_names = classnames(G;TeX=true)
    weyl_classes_names = fromTeX.(Ref(rio()),weyl_classes_TeX_names)

    # Remove rational orbit labels
    rational_label_indices = [ x for g in groups for x in g[2:end] ] # Find labels which need to be removed
    keep = setdiff(1:length(rational_unipotent_classes_names), rational_label_indices)
    geometric_unipotent_classes_names = rational_unipotent_classes_names[keep] 

    # Formatting table
    repr_data = xrepr.(Ref(rio()),CycPol.(intersection_numbers_scaled_summed))
    println("")
    println("Rows labelled by conj classes C∈ConjCl(W)")
    println("Columns labelled by unipotent classes C⊆G(Fq)")
    println("")
    showtable(repr_data;col_labels=geometric_unipotent_classes_names,row_labels=weyl_classes_names,rows_label="|BwB∩C|")
end