# Eulerian distance for internal use
BindGlobal("__cryst__EuclideanDistance",
    function(x, y)
        return Sqrt( Float(x[1]-y[1])^2 + Float(x[2]-y[2])^2 + Float(x[3]-y[3])^2);
    end
);

# Eulerian distance for internal use
BindGlobal("__cryst__EuclideanNorm",
    function(x)
        local res, minDist, i;
        res:=0;
        for i in [1..Length(x)] do
            res:=res+x[i]^2;
        od;
        return Sqrt(res);
    end
);


InstallGlobalFunction( DirichletCellForFiniteWord,
    function( vector, length, generatingSet)
        local i, j, k, elementsOfLenghtL, gen, elementsInOrbit, e, hyperplanes, y, h, l, el, fundDomSurface, vif, coords, polymakeObj, vifTriangulated, facett, nextEdge, edges, facet, orderedFacets, newFacet, currVertex, possVertex, remainingVertices, oldElems, triangulatedFacets, pr;

        elementsOfLenghtL:=DuplicateFreeList(generatingSet);

        # check that vector has trailing 1 and if not append it
        if Length(vector) in [3,4] then
            if Length(vector) = 4 and not vector[4] = 1 then
                Error("last entry of vector to construct dirichlet cell around needs to end with trailing 1 or have three entries");
            fi; 
            if Length(vector) = 3 then
                Add(vector, 1);
            fi;
        else
            Error("vector needs to either have three entries or four entries including a trailing 1");
        fi; 

        for i in [1..Length(vector)] do
            vector[i] := Rat(vector[i]);
        od;

        # generate all elements induced by words of maximum length length
        for i in [1..length-1] do
            for gen in generatingSet do
                oldElems:=StructuralCopy(elementsOfLenghtL);
                for el in oldElems do
                    if not el*gen in elementsOfLenghtL then
                        Add(elementsOfLenghtL, el*gen);
                        #Print("added element ", el*gen, "\n");
                    fi;
                    if not el*(gen^(-1)) in elementsOfLenghtL then
                        Add(elementsOfLenghtL, el*(gen^(-1)));
                        #Print("added element ", el*(gen^(-1)), "\n");
                    fi;
                od;
            od;
        od;
        #Print("elements generated\n");

        # generate all points of x^e for e \in elementsOfLenghtL
        elementsInOrbit := [];
        for e in elementsOfLenghtL do
            Add(elementsInOrbit, (vector*e));
        od;

        #Print("Orbits generated\n");

        # calculate the hyperplanes necessary for the dirichlet construction
        # we cut of the last coordinate as it is always 0 because of the trailing 1 from the isometry operations
        # using the formulas presented here: https://math.stackexchange.com/questions/2858815/understanding-formula-for-hyperplanes
        hyperplanes:=[];
        for y in elementsInOrbit do
            h:=[];
            h[1]:=1/2*(__cryst__EuclideanNorm((y))^2 - __cryst__EuclideanNorm(Rat(vector))^2);
            for i in [1..Length(vector)-1] do
                h[i+1]:=Rat(y[i]-vector[i]);
            od;
            Add(hyperplanes, h);
        od;

        #Print("hyperplanes generated\n");

        # remove zeroes, as they are not allowed
        while not Position(hyperplanes, [0.,0.,0.,0.]) = fail do
            Remove(hyperplanes, Position(hyperplanes, [0.,0.,0.,0.]));
        od;

        # Print(hyperplanes);

        # get the vertex coordinates with polymake
        polymakeObj:=CreatePolymakeObject();;
        AppendInequalitiesToPolymakeObject(polymakeObj, hyperplanes);;
        coords:=Polymake(polymakeObj, "VERTICES");;
        # get vertices in faces as well
        vif:=Polymake(polymakeObj, "VERTICES_IN_FACETS");;

        # vertices in faces are not necessary triangulated
        vifTriangulated:=[];
        edges:=[];
        for facet in vif do
            for i in [1..Size(facet)] do
                for j in [i..Size(facet)] do
                    if not i=j then
                        Add(edges, [facet[i], facet[j]]);
                    fi;
                od;
            od;
        od;

        # order the facets
        orderedFacets:=[];
        for facet in vif do
            if Size(facet) = 3 then
                Add(orderedFacets, facet);
            else
                newFacet:=[facet[1]];
                remainingVertices:=StructuralCopy(facet);
                Remove(remainingVertices, 1);
                while not IsEmpty(remainingVertices) do
                    currVertex:=Last(newFacet);
                    nextEdge:=[];
                    for possVertex in remainingVertices do
                        if Size(Filtered(edges, e -> e=[possVertex, currVertex] or e=[currVertex, possVertex]))>1 then
                            nextEdge:=Filtered(edges, e -> e=[possVertex, currVertex] or e=[currVertex, possVertex])[1];
                        fi;
                    od;
                    if nextEdge = [] then
                        Error("could not find new edge, currently trying from vertex ", currVertex, "\n");
                    fi;
                    if nextEdge[1] = currVertex then
                        Add(newFacet, nextEdge[2]);
                        Remove(remainingVertices, Position(remainingVertices, nextEdge[2]));
                    else
                        Add(newFacet, nextEdge[1]);
                        Remove(remainingVertices, Position(remainingVertices, nextEdge[1]));
                    fi;
                od;
                Add(orderedFacets, newFacet);
            fi;
        od;

        # triangulate the ordered faces
        triangulatedFacets:=[];
        for facet in orderedFacets do
            if not Size(facet) = 3 then
                for i in [2..(Size(facet)-1)] do
                    Add(triangulatedFacets, [facet[1], facet[i], facet[i+1]]);
                od;
            else
                Add(triangulatedFacets, facet);
            fi;
        od;

        # plot with GAPic
        # fundDomSurface:=SimplicialSurfaceByVerticesInFaces(triangulatedFacets);
        # pr:=SetVertexCoordinates3D(fundDomSurface, coords);
        # DrawComplexToJavaScript(fundDomSurface, Concatenation("dirichlet-cell-", String(length)), pr);

        Print("dirichlet cell has volume of: ", Float(Polymake(polymakeObj, "VOLUME")), "\n\n");

        # return [fundDomSurface, triangulatedFacets, coords];
        return [triangulatedFacets, coords];
        # return 0;
    end 
);