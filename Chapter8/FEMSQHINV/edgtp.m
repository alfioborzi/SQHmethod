function [type] = edgtp(EDGE,element,vertex);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    get edge type from list of edges.
%
% notes: if edge is not in list, then it is an interior edge, 
%        in which case we return 0.
%
% output: type=0 ==> edge not in EDGE array (interior edge)
%         type=1 ==> dirichlet boundary edge
%         type=2 ==> neumann boundary edge
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%% default return
   type = 0;

   %%% make sure we have SOME edges...
   if (EDGE == [])
      fprintf('edgetype error: edge array is null...\n');
      return;
   end;

   %%% consider only edges belonging to the specified element
   matching_elements = find(EDGE(:,1)==element);
   [size1,size2] = size(matching_elements);
   if ( size1 == 0 ) | ( size2 == 0 )
      return;
   end;
   EDGE_smaller = EDGE(matching_elements,:);

   %%% now see which of the remaining edges start with the specified vertex
   vertex_list = EDGE_smaller(:,2);
   matches = find(vertex_list==vertex);
   [size1,size2] = size(matches);
   if ( size1 == 1) & ( size2 == 1 )
      type = EDGE_smaller(matches,3);
      return;
   end;
   if ( size1 > 1 ) | ( size2 > 1)
      fprintf('edgetype error: edge occurs twice in edge array...\n');
      return;
   end;

