FUNCTION CONVEXHULL_volume, x, y, z, verts = verts, conns = conns

; PURPOSE:
;   Calculate the volume of convex hull structure for a 3D-structure
;   More details about mesh in IDL refer to book "Modern IDL: A Guide
;   to IDL Programming" by Michael Galloy


; INPUT:
;  x, y, z - input x, y, z array of a 3D-structure

; 0.Construct the vertex of each pixel
;   (1) build vertexs
  npixels = n_elements(x)
  xver = make_array(npixels*8, /float)
  yver = make_array(npixels*8, /float)
  zver = make_array(npixels*8, /float)
  for i=0, npixels-1 do begin
	  ; celling vertex
      xver[i*8]=x[i]
      xver[i*8+1]=x[i]+1.0
      xver[i*8+2]=x[i]+1.0
      xver[i*8+3]=x[i]

      yver[i*8]=y[i]
      yver[i*8+1]=y[i]
      yver[i*8+2]=y[i]+1.0
      yver[i*8+3]=y[i]+1.0

      zver[i*8:i*8+3] = z[i]
	  ; bottome vertex z = 0
      xver[i*8+4]=x[i]
      xver[i*8+5]=x[i]+1.0
      xver[i*8+6]=x[i]+1.0
      xver[i*8+7]=x[i]

      yver[i*8+4]=y[i]
      yver[i*8+5]=y[i]
      yver[i*8+6]=y[i]+1.0
      yver[i*8+7]=y[i]+1.0

      zver[i*8+4:i*8+7] = 0.0

  endfor


;   (2) remove duplicates
  xsize = max(x) - min(x) + 2
  ysize = max(y) - min(y) + 2
  zsize = max(z) - min(z) + 2
  l3d = xver + yver*xsize + zver*xsize*ysize
  uniq_ind = uniq(l3d, sort(l3d))

  xver = xver[uniq_ind]
  yver = yver[uniq_ind]
  zver = zver[uniq_ind]


; 1.Construct the Delaunay triangulation
;   and the Voronoi diagram.
  QHULL, xver, yver,zver, triangle, BOUNDS=bounds

; 2. Calculate the verts and connectivity array
;   2.1 verts
  verts = [transpose(xver), transpose(yver), transpose(zver)]

;   2.2 connectivity array, more information refer to 
;   https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-polygon-mesh
;   More information about connectivity - conns refer to below website:
;   file:////Applications/exelis/idl85/help/online_help/IDL/idl.htm#Object%20Classes/Graphics/IDLgrPolygon.htm?Highlight=polygons
  conn = make_array(n_elements(triangle[0,*]), /INTEGER) + 3
  conns = [transpose(conn), triangle]
  conns = reform(conns, n_elements(conns))

  if MESH_ISSOLID(conns) eq 1 then $
     volume = MESH_VOLUME(verts, conns) else $
     volume = -1

 ; 3. Plot convex hull structure 
 ;window, xsize = 300, ysize=300, /free
 ;scale3, xrange=[0,100], yrange=[0,100], zrange=[0,100]
 ;p = polyshade(verts, conns, /t3d)
 ;tv, p


  return, volume


END
