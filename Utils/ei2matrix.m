function [EIm, EImView] = ei2matrix(oneEI)
% Converts EI data sorted by electrode from 512-electrode array to a 2d matrix
%      input:   oneEI:   ei data of size = [1,512]
%
%    outputs:   EIm:     ei data reshaped into a matrix of size [16,32],
%                        corresponding roughly to electrode locations 
%                        (ignores hexagonal arrangement) 
%               EImView: ei data reshaped to a matrix of size [32,65],
%                        corresponding to electrode locations and including
%                        hexagonal electrode arrangement, for displaying ei 
%                        using imagesc(). 4 pixels per electrode
%                       
% L Grosberg Jan 2014

% Load electrode positions as matrix 
fname = fullfile(fileparts(mfilename('fullpath')),'../resources/array_matrix_id510'); 
h = load(fname); 
electrodeMatrix = h.array_matrix_id510; clear h; 

EIm = oneEI(electrodeMatrix(:)); 
EIm = reshape(EIm, size(electrodeMatrix)); %figure; imagesc(EIm);axis image; 

% Reshape the matrix so it can be viewed accurately as an image
EImView = upsample(EIm',2)'; 
EImView(:,2:2:end) = EImView(:,1:2:end); 
% Shift even rows by one
EImView = cat(2,EImView,zeros(size(EImView,1),1)); %Pad all rows with zero
for odds = 1:2:size(EImView,1); 
    EImView(odds,:) = circshift(EImView(odds,:),[0 1]); 
end
EImView            = upsample(EImView,2); 
EImView(2:2:end,:) = EImView(1:2:end,:); 

%figure; imagesc(EIm,[0 max(EIm(:))/2]); axis image; colorbar; title('EI'); 
end