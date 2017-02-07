% Function to make an x or a square for plotting on the image
%
% Input:
% sizeOfShape The dimension of the shape
% shape The shape you want to make 'square' or 'x'
% thickness The thickness of the plot 'thick' or 'thin'
%
% Output:
% shapeMatrix The matrix representing the square or the x
%
% John Ginn
% Created: 7/22/16
% Modified: 7/22/16

function [shapeMatrix] = makeShape(sizeOfShape,shape,thickness)
shapeMatrix = zeros(sizeOfShape,sizeOfShape);

if strcmp(shape,'x') == 1
    % make the 'x'
    for step = 1:sizeOfShape
        shapeMatrix(step,step) = 1;
        shapeMatrix(step,sizeOfShape - step + 1) = 1;
        % make the shape thicker
        if strcmp(thickness,'thick') == 1
            % stay in bounds for extra thickness of shape
            if((step + 1) <= sizeOfShape)
                shapeMatrix(step,(step + 1)) = 1;
            end
            if((step - 1) >= 1)
                shapeMatrix(step,(step - 1)) = 1;
            end
            if((sizeOfShape - step + 2) <= sizeOfShape)
                shapeMatrix(step,sizeOfShape - step + 2) = 1;
            end
            if((sizeOfShape - step) >= 1)
                shapeMatrix(step,sizeOfShape - step) = 1;
            end
        end
    end
elseif strcmp(shape,'square')
    % make the square
    for step = 1:sizeOfShape
        if (step == 1)||(step == sizeOfShape);
            % make the top row or bottom row
           for stepCol = 1:sizeOfShape
               shapeMatrix(step,stepCol) = 1;
           end
        else
            % add a point on either size
        shapeMatrix(step,1) = 1;
        shapeMatrix(step,sizeOfShape) = 1;
        end
        
%        % make the shape thicker
%        if strcmp(thickness,'thick') == 1        
%         % make square thicker
%         if (step == 2)||(step == sizeOfShape-1);
%             % make the top row or bottom row
%             for stepCol = 1:sizeOfShape
%                 shapeMatrix(step,stepCol) = 1;
%             end
%         else
%             % add a point on either size
%             shapeMatrix(step,2) = 1;
%             shapeMatrix(step,sizeOfShape-1) = 1;
%         end
%        end
    end
end

end
