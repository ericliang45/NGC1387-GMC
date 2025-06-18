PRO qhull_example2
 
   ; Create a cube vertix
   x=[0, 1, 1, 0, 0, 1, 1, 0, 2, 2, 1, 2, 2, 1]
   y=[0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1]
   z=[0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 2, 2, 2, 2]
 
   ; Construct the Delaunay triangulation
   ; and the Voronoi diagram.
   ;QHULL, x, y, triangle, /DELAUNAY, $
   QHULL, x, y, z, triangle, BOUNDS=bounds
;      /DELAUNAY, CONNECTIVITY=connectivity,$
;      VDIAGRAM=vdiagram, VVERTICES=vvert, VNORM=vnorm
 
    print,'triangle',triangle
    print, 'bounds',bounds

	verts = [transpose(x), transpose(y), transpose(z)]


    conn1 = make_array(n_elements(triangle[0,*]), /INTEGER) + 3
    ;conn2 = make_array(n_elements(triangle[0,*]), /INTEGER) + 0
    conns = [transpose(conn1), triangle]
    conns = reform(conns, n_elements(conns))

	print,'conns',conns
  
    print, 'MESH_ISSOLID(conns)',MESH_ISSOLID(conns)
    print, 'MESH_volume(verts, conns)',MESH_volume(verts, conns)

	x1 = [0, 1]
	y1 = [0, 0]
	z1 = [1, 2]
	print, 'convexhull_volume(x1,y1,z1)',convexhull_volume(x1,y1,z1)
 
 
END
