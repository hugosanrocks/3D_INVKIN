       program test_delaunay_mod

       !**************************************************************!
       !
       !! MAIN is the main program for TABLE_DELAUNAY.
       !
       !  Discussion:
       !
       !    TABLE_DELAUNAY computes the Delaunay triangulation of a TABLE dataset.
       !
       !    The dataset is simply a set of points in the plane.
       !
       !    Thus, given a set of points V1, V2, ..., VN, we apply a standard 
       !    Delaunay triangulation.  The Delaunay triangulation is an organization 
       !    of the data into triples, forming a triangulation of the data, with
       !    the property that the circumcircle of each triangle never contains
       !    another data point.  
       !
       !  Usage:
       !
       !    table_delaunay prefix
       !
       !    where:
       !
       !    * prefix_nodes.txt,     the node coordinates (input).
       !    * prefix_elements.txt,  the nodes that make up each triangle (output).
       !
       !  Licensing:
       !
       !  Modified:
       !
       !    12 March 2018
       !
       !  Author:
       !
       !    Hugo Sanchez
       !
       !  Usage:
       !
       !    call dtris2 ( node_num, node_xy, triangle_num, triangle_node, &
       !                  triangle_neighbor )
       !    input:  node_num      	total number of nodes
       !            node_xy       	node coordinates
       !    output: triangle_num    	total number of triangles
       !            triangle_node   	node number of 3 vertex of triangles
       !            triangle_neighbor 	triangle number of neighbor
       !
         use delaunay
         implicit none
         integer ( kind = 4 ) node_dim
         integer ( kind = 4 ) node_num
         real ( kind = 4 ), allocatable, dimension ( :, : ) :: node_xy
         integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
         integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node
         integer ( kind = 4 ) triangle_num
         integer ( kind = 4 ) triangle_order
         integer ( kind = 4 ) iunit, i
       
         node_dim = 2
         node_num = 5745
       
         allocate ( node_xy(1:node_dim,1:node_num) )
       
         !Read node coordinates
         call read_nodes(node_num,node_xy)
 
         !Determine the Delaunay triangulation.
         triangle_order = 3
       
         allocate ( triangle_node(triangle_order,3*node_num) )
         allocate ( triangle_neighbor(triangle_order,3*node_num) )
       
         call dtris2 ( node_num, node_xy, triangle_num, triangle_node, &
                       triangle_neighbor )
       
         write ( *, '(a,i8)' ) '  Number of triangles is ', triangle_num
       
         !Write the triangulation to a file.
         call write_triangles(triangle_num,triangle_node)

         !Free memory.
         deallocate ( node_xy )
         deallocate ( triangle_node )
         deallocate ( triangle_neighbor )
       
         stop
       endprogram test_delaunay_mod
