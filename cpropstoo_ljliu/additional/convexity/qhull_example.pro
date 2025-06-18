PRO qhull_example
 
   ; Create a collection of random points.
   n = 20
   seed = 15
   x = RANDOMU(seed, n)
   y = RANDOMU(seed, n)
 
   ; Construct the Delaunay triangulation
   ; and the Voronoi diagram.
   ;QHULL, x, y, triangle, /DELAUNAY, $
   QHULL, x, y, triangle, BOUNDS=bounds
;      /DELAUNAY, CONNECTIVITY=connectivity,$
;      VDIAGRAM=vdiagram, VVERTICES=vvert, VNORM=vnorm
 
   ; Plot our input points.
   PLOT, [-0.1, 1.1], [-0.1, 1.1], /NODATA, $
      XSTYLE=4, YSTYLE=4
   PLOTS, x, y, PSYM=4
   for i = 0, 19 do begin
     xyouts,x[i],y[i],strtrim(i,2), charsize=2.
   endfor
 
   ; Plot the Voronoi diagram.
;   FOR i=0,N_ELEMENTS(vdiagram[2,*])-1 DO BEGIN
;      vdiag = vdiagram[*, i]
;      j = vdiag[2] + 1
;      ; Bounded or unbounded?
;      IF (j gt 0) THEN BEGIN   ; Bounded.
;         PLOTS, vvert[*, vdiag[2:3]], PSYM=-5
;         CONTINUE
;      ENDIF
; 
;      ; Unbounded, retrieve starting vertex.
;      xystart = vvert[*, vdiag[3]]
;      ; Determine the line equation.
;      ; Vnorm[0]*x + Vnorm[1]*y + Vnorm[2] = 0
;      slope = -vnorm[0,-j]/vnorm[1,-j]
;      intercept = -vnorm[2,-j]/vnorm[1,-j]
; 
;      ; Need to determine the line direction.
;      ; Pick a point on one side along the line.
;      xunbound = xystart[0] + 5
;      yunbound = slope*xunbound + intercept
; 
;      ; Find the closest original vertex.
;      void = MIN( (x-xunbound)^2 + (y-yunbound)^2, idx)
;      ; By definition of Voronoi diagram, the line should
;      ; be closest to one of the bisected points. If not,
;      ; our line went in the wrong direction.
;      IF (idx ne vdiag[0] && idx ne vdiag[1]) THEN BEGIN
;         xunbound = xystart[0] - 5
;         yunbound = slope*xunbound + intercept
;      ENDIF
; 
;      PLOTS, [[xystart], [xunbound, yunbound]]
; 
;   ENDFOR

    for i=0, n_elements(triangle[0,*])-1 do begin
        plots, [x[triangle[*,i]],x[triangle[0,i]]], $
           [y[triangle[*,i]],y[triangle[0,i]]], psym=-5
    endfor

    print,'triangle',triangle
    print, 'bounds',bounds


    conn1 = make_array(n_elements(triangle[0,*]), /INTEGER) + 4
    conn2 = make_array(n_elements(triangle[0,*]), /INTEGER) + 3
    conns = [transpose(conn2), transpose(conn1), triangle]
    conns = reform(conns, n_elements(conns))

	print,'conns',conns
  
    print, 'MESH_ISSOLID(conns)',MESH_ISSOLID(conns)
 
 
END
