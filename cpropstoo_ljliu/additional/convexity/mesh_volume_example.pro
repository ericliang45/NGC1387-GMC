pro mesh_volume_example
    n = 40
    mid = n / 2
    radius = n / 3.5
    vol = BYTARR(n,n,n)
    for i=0, n-1 do begin
        for j=0, n-1 do begin
            for k=0, n-1 do begin
                vol[i,j,k] = SQRT( (i-mid)^2 + (j-mid)^2 + (k-mid)^2 ) + 0.5
            endfor
        endfor
    endfor

    ISOSURFACE, vol, radius, v, c

	print,'size(v)',size(v)

    print, "Ideal volume", 4./3.*!pi*radius^3
    print, "MESH_ISSOLID(c)",MESH_ISSOLID(c)
	print,'c',c
    print, "Isosurface volume", MESH_VOLUME(v, c)

    INTERVAL_VOLUME, vol, 0, radius, tet_verts, tet_conn
    print, "Tetrahedral volume", TETRA_VOLUME(tet_verts, tet_conn)
    surf_conn = TETRA_SURFACE(tet_verts, tet_conn)
    print, MESH_ISSOLID(surf_conn)
    print, "Tetrahedral volume by surface", MESH_VOLUME(tet_verts, surf_conn)
end
