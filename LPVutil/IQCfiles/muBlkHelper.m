function blk = muBlkHelper(G)
% Returns NBLK-by-2 array representing the uncertainty "block" structure of
% G represented in the simple form used by MUSSV.
%
% Related to private function, blkstruct2N2, which is called within MUSSV,
% and allows a user to call MUSSV using the 3rd argument from LFTDATA as
% the block structure.
[~,~,a3,~] = lftdata(G);
nblk = numel(a3);
blk = zeros(nblk,2);
for i=1:nblk
   switch a3(i).Type
      case 'ureal'
         blk(i,:) = [-a3(i).Occurrences 0];
      case {'ultidyn', 'ucomplexm'}
         szB = a3(i).Size;
         if a3(i).Occurrences==1
            blk(i,:) = szB;
         elseif all(szB==1)
            blk(i,:) = [a3(i).Occurrences 0];
         else
            error('Repeated, non-scalar complex or ultidyn blocks cannot be represented');
         end
      case 'ucomplex';
         blk(i,:) = [a3(i).Occurrences 0];
      case 'udyn'
         error('UDYN blocks cannot be represented');
      otherwise
   end
end
   