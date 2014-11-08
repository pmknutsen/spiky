function [counts,x_inds,y_inds] = hist2d(data_mtx, varargin)
% HIST2D  2-Dimensional Density Histogram.
% 
% COUNTS = HIST2D(DATA_MTX, BINS), where 'data_mtx' is an [m x n] input
%    matrix produces a 'counts' matrix ([bins x n]) where each column is
%    a histogram of the corresponding column of 'data_mtx'.  If the number
%    of bins is not specified in the 'bins' argument, it defaults to 100.
%
% [COUNTS,X_INDS,Y_INDS] = hist2d(DATA_MTX) also returns the bin centers
%    and column indices as vectors, so the histogram can be visualized with
%    IMAGESC(X_INDS,Y_INDS,COUNTS).
%
% HIST2D(...,'log') computes the log of the counts (both for display and
%    in the returned 'counts').
% 
% A special case arises when 'data_mtx' is an [m x 2] matrix.  In this
%    case, the 'counts' output matrix will be of size [bins x bins] and
%    represents the densities of points a discretized scatter plot of the
%    first column vs. the second.
%
% HIST2D(...) without output arguments produces an image of the density.

bins = 100;
logarg = 0;
for arg = 1:length(varargin)
    if (strcmp(class(varargin{arg}),'char'))
        if (~strcmp(lower(varargin(arg)),'log'))
            error('Unknown command option.');
        else
            logarg = 1;
        end
    elseif (strcmp(class(varargin{arg}), 'double') & (numel(varargin{arg}) == 1))
        bins = round(varargin{arg});
    else
        error('Unknown argument');
    end
end

dims = size(data_mtx);

if (dims(2) < 2)
    error('Data matrix must be [m x n] with n greater than 1.');
elseif (dims(2) == 2)  % scatter plot density
    % If we're scatter plotting, first separately scale each column
    [data1,min1,max1] = rescale(data_mtx(:,1),1,bins);
    [data2,min2,max2] = rescale(data_mtx(:,2),1,bins);
    
    % The bin centers are built based on range of the unscaled data.
    x_inds = linspace(min1,max1,bins);
    y_inds = linspace(min2,max2,bins);
    
    % We do the histogramming by using the (scaled & rounded) values of
    % the data columns as row and column indices and relying on 'sparse'
    % to keep count for us (the use of the sparse matrix is incidental;
    % we do this just because its a builtin way of counting).
    counts = full(sparse(round(data2),round(data1),1,bins,bins));
    
else % data density
    % Here all of the data is assumed to come from the same distribution so
    % it is all rescaled at once.
    [data_mtx,oldmin,oldmax] = rescale(data_mtx,1,bins);
    data_mtx = round(data_mtx);
    
    % The row bin centers in again taken from the data but the column
    % indices are just the data column indices.
    y_inds = linspace(oldmin,oldmax,bins);
    x_inds = 1:dims(2); 
    
    % Again, we use the 'sparse' function to histogram for us, but here
    % we use the data value as the row index while preserving the column
    % index.  This gives us a column by column histogram.
    colvector  = reshape(repmat(x_inds, [dims(1),1]), [numel(data_mtx),1]);
    datavector = round(reshape(data_mtx, [numel(data_mtx),1]));
    counts = full(sparse(double(datavector),colvector,1,bins,dims(2)));
end

if (logarg)
    old = warning('off');
    counts = log(counts);
    warning(old);
end

%% Plot the density if no outputs requested
if (nargout == 0)
%     h = surf(x_inds,y_inds,counts);       % these lines draw a surface instead of an image ...
%     set(h, 'AmbientStrength', 0.6); view([45,80]); nice;
     imagesc(x_inds, y_inds, counts); axis xy;
     clear counts x_inds y_inds  % clear these so nothing is dumped to output
end
