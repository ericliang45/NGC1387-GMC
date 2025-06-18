FUNCTION CONVEXHULL_AREA, x, y

; PURPOSE:
;   Calculate the area of convex hull structure for a 2D-structure

; INPUT:
;  x, y - input x, y array of a 2D-structure

; 0.Construct the vertex of each pixel
;   (1) build vertexs
  npixels = n_elements(x)
  xver = make_array(npixels*4, /float)
  yver = make_array(npixels*4, /float)
  for i=0, npixels-1 do begin
	  xver[i*4]=x[i]
	  xver[i*4+1]=x[i]+1.0
	  xver[i*4+2]=x[i]+1.0
	  xver[i*4+3]=x[i]

	  yver[i*4]=y[i]
	  yver[i*4+1]=y[i]
	  yver[i*4+2]=y[i]+1.0
	  yver[i*4+3]=y[i]+1.0
  endfor

;   (2) remove duplicates
  xsize = max(x) - min(x) + 2
  ysize = max(y) - min(y) + 2
  l2d = xver + yver*xsize 
  uniq_ind = uniq(l2d)


  xver = xver[uniq_ind]
  yver = yver[uniq_ind]

; 1.Construct the Delaunay triangulation
;   and the Voronoi diagram.
  QHULL, xver, yver, triangle, BOUNDS=bounds, $
	        VDIAGRAM=vdiagram, VVERTICES=vvert, VNORM=vnorm

; 2. Calculate the area of the convex hull structure
   area = 0.
   for i=0, n_elements(triangle[0,*])-1 do begin
      area = area + poly_area(xver[triangle[*,i]],yver[triangle[*,i]])
   endfor

   return, area


END
