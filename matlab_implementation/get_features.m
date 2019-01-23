function x = get_features(im, features, cell_size, cos_window)
%GET_FEATURES
%   Extracts dense features from image.
%
%   X = GET_FEATURES(IM, FEATURES, CELL_SIZE)
%   Extracts features specified in struct FEATURES, from image IM. The
%   features should be densely sampled, in cells or intervals of CELL_SIZE.
%   The output has size [height in cells, width in cells, features].
%
%   To specify HOG features, set field 'hog' to true, and
%   'hog_orientations' to the number of bins.
%
%   To experiment with other features simply add them to this function
%   and include any needed parameters in the FEATURES struct. To allow
%   combinations of features, stack them with x = cat(3, x, new_feat).
%
%   Joao F. Henriques, 2014
%   http://www.isr.uc.pt/~henriques/


	if features.hog,
		%HOG features, from Piotr's Toolbox
		x = double(fhog(single(im) / 255, cell_size, features.hog_orientations));
		x(:,:,end) = [];  %remove all-zeros channel ("truncation feature")
	end
	
	if features.gray,
		%gray-level (scalar feature)
		x = double(im) / 255;
		
		x = x - mean(x(:));
	end
	
	%process with cosine window if needed
    
    if isfield(features, 'deep') && features.deep
        x = impreprocess(single(im));
%         caffe('reshape_input', 'solver', [0, 1, size(x, 3), size(x, 2), size(x, 1)]);
%         x = caffe('forward', {x});
        caffe('set_input_dim', 'DNNL', [0, 1, size(x, 3), size(x, 2), size(x, 1)]);
        x = caffe('forward', 'DNNL', {x});
        x = permute(x{1}, [2, 1, 3]);
        x = x / 1e3;
    end
    
	if ~isempty(cos_window),
		x = bsxfun(@times, x, cos_window);
	end
	
end
