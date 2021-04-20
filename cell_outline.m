function outline = cell_outline(Im)
%CELL_OUTLINE - returns outline of the sell
%   Im - cell matrix (1s and 0s)

outline = zeros(size(Im)+2);
Im = [zeros(1,size(Im,1)+2);zeros(size(Im,1),1), Im, zeros(size(Im,1),1);zeros(1,size(Im,1)+2)];
outline(2:end-1,2:end-1) = Im(1:end-2,2:end-1) + Im(3:end,2:end-1) + ...
    Im(2:end-1,1:end-2) + Im(2:end-1,3:end);
outline = outline - 4*Im;
o = outline>0;
o = o(2:end-1, 2:end-1);
outline = o;

end

