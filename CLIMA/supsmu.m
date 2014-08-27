function smo = supsmu(x,y,varargin)
%supsmu: Smoothing of scatterplots using Friedman's supersmoother algorithm.
%
% Syntax:
%   Y_SMOOTH = supsmu(X,Y,'PropertyName',PropertyValue,...)
%
% Inputs:
%   X, Y are same-length vectors.
%
% Output
%   Y_SMOOTH is a smoothed version of Y.
%
% Example,
%   x = linspace(0,1,201);
%   y = sin(2.5*x) + 0.05*randn(1,201);
%   smo = supsmu(x,y);
%   plot(x,y,'o',x,smo)
%
% The supersmoother algorithm computes three separate smooth curves from
% the input data with symmetric spans of 0.05*n, 0.2*n and 0.5*n, where n
% is the number of data points.  The best of the three smooth curves is
% chosen for each predicted point using leave-one-out cross validation. The
% best spans are then smoothed by a fixed-span smoother (span = 0.2*n) and
% the prediction is computed by linearly interpolating between the three
% smooth curves.  This final smooth curve is then smmothed again with a
% fixed-span smoother (span = 0.05*n).
%
% According to comments by Friedman, "For small samples (n < 40) or if
% there are substantial serial correlations between observations close in
% x-value, then a prespecified fixed span smoother (span > 0) should be
% used.  Reasonable span values are 0.2 to 0.4."
%
%
% The following optional property/value pairs can be specified as arguments
% to control the indicated behavior:
%
%   Property    Value
%   ----------  ----------------------------------------------------------
%   Weights     Vector of relative weights of each data point.  Default is
%               for all points to be weighted equally.
%
%   Span        Sets the width of a fixed-width smoothing operation
%               relative to the number of data points, 0 < SPAN < 1.
%               Setting this to be non-zero disables the supersmoother
%               algorithm.  Default is 0 (use supersmoother).
%
%   Period      Sets the period of periodic data.  Default is Inf
%               (infinity) which implies that the data is not periodic.
%               Can also be set to zero for the same effect.
%
%   Alpha       Sets a small-span penalty to produce a greater smoothing
%               effect.  0 < Alpha < 10, where 0 does nothing and 10
%               produces the maximum effect.  Default = 0.
%
%   Unsorted    If the data points are not already sorted in order of the X
%               values then setting this to true will sort them for you.
%               Default = false.
%
% All properties names are case-insensitive and need only be unambiguous.
% For example,
%
%   Y_SMOOTH = supsmu(X,Y,'weights',ones(n,1),'per',2*pi)
%
% is valid usage.

% Friedman, J. H. (1984). A Variable Span Smoother. Tech. Rep. No. 5,
% Laboratory for Computational Statistics, Dept. of Statistics, Stanford
% Univ., California.

% Version: 1.0, 12 December 2007
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Input checks.
error(nargchk(2,Inf,nargin))

% x and y must be vectors with same number of points (at least 5).
sy = size(y);
n = length(y);
if ~isvector(x) || ~isvector(y) || length(x) ~= n || n < 5
	error('X and Y must be equal-length vectors of at least 5 points.')
end

% Define properties and set default values.
prop.weights = [];
prop.span = 0;
prop.period = Inf;
prop.alpha = 0;
prop.unsorted = false;

% Process inputs and set prop fields.
properties = fieldnames(prop);
arg_index = 1;
while arg_index <= length(varargin)
	arg = varargin{arg_index};
	if ischar(arg)
		prop_index = find(strncmpi(arg,properties,length(arg)));
		if length(prop_index) == 1
			prop.(properties{prop_index}) = varargin{arg_index + 1};
		else
			error('Property ''%s'' does not exist or is ambiguous.',arg)
		end
		arg_index = arg_index + 2;
	elseif isstruct(arg)
		arg_fn = fieldnames(arg);
		for i = 1:length(arg_fn)
			prop_index = find(strncmpi(arg_fn{i},properties,...
				length(arg_fn{i})));
			if length(prop_index) == 1
				prop.(properties{prop_index}) = arg.(arg_fn{i});
			else
				error('Property ''%s'' does not exist or is ambiguous.',...
					arg_fn{i})
			end
		end
		arg_index = arg_index + 1;
	else
		error(['Properties must be specified by property/value pairs',...
			' or structures.'])
	end
end

% Validate Weights property.
if isempty(prop.weights)
elseif length(prop.weights) == 1
	prop.weights = [];
elseif isvector(prop.weights) && length(prop.weights) == n
	prop.weights = prop.weights(:);
else
	error('Weights property must be a vector of the same length as X and Y.')
end

% Validate Span property.
if ~isscalar(prop.span) || prop.span < 0 || prop.span >= 1
	error('Span property must be a scalar, 0 <= span < 1.')
end

% Validate Periodic property.
if ~isscalar(prop.period) || prop.period < 0
	error('Periodic property must be a scalar >= 0.')
end
if isinf(prop.period)
	prop.period = 0;
end

% Validate Alpha property.
if isscalar(prop.alpha)
	prop.alpha = min(max(prop.alpha,0),10);
else
	error('Alpha property must be a scalar.')
end

% Validate Unsorted property.
if ~isscalar(prop.unsorted)
	error('Unsorted property must be a scalar.')
end


% Select one of four smooth functions.  Each smooth function has been
% speed optimized for the specific conditions.
smooth_fcn_selector = 2*isempty(prop.weights) + (prop.period == 0);
switch smooth_fcn_selector
	case 0 % use weights, periodic
		smooth = @smooth_wt_per;
	case 1 % use weights, aperiodic
		smooth = @smooth_wt_aper;
	case 2 % no weights, periodic
		smooth = @smooth_per;
	case 3 % no weights, aperiodic
		smooth = @smooth_aper;
end

% Make x and y into column vectors and sort if necessary.
x = x(:);
y = y(:);
if prop.unsorted
	[x,order] = sort(x);
	y = y(order);
	if ~isempty(prop.weights)
		prop.weights = prop.weights(order);
	end
end

% If prop.span > 0 then we have a fixed span smooth.
if prop.span > 0
	smo = smooth(x,y,prop.weights,prop.span,prop.period);
	smo = reshape(smo,sy);
	if prop.unsorted
		smo(order) = smo;
	end
	return
end

spans = [0.05;0.2;0.5];
nspans = length(spans);

% Compute three smooth curves.
smo_n = zeros(n,nspans);
acvr_smo = zeros(n,nspans);
for i = 1:nspans
	[smo_n(:,i),abs_cv_res] = smooth(x,y,prop.weights,spans(i),prop.period);
	acvr_smo(:,i) = smooth(x,abs_cv_res,prop.weights,spans(2),prop.period);
end

% Select which smooth curve has smallest error using cross validation.
[resmin,index] = min(acvr_smo,[],2);
span_cv = spans(index);

% Apply alpha.
if prop.alpha ~= 0
	small = 1e-7;
	tf = resmin < acvr_smo(:,3) & resmin > 0;
	span_cv(tf) = span_cv(tf) + (spans(3) - span_cv(tf)).* ...
		max(small,resmin(tf)./acvr_smo(tf,3)).^(10 - prop.alpha);
end

% Smooth span_cv and clip at spans(1) and spans(end).
smo_span = smooth(x,span_cv,prop.weights,spans(2),prop.period);
smo_span = max(min(smo_span,spans(end)),spans(1));

% Interpolate each point.
% The block of code below does the same thing as this, but much faster:
% smo_raw = zeros(n,1);
% for i = 1:n
% 	smo_raw(i) = interp1(spans,smo_n(i,:),smo_span(i));
% end
try
	bin = sum(bsxfun(@ge,smo_span,spans(1:end-1)'),2);
catch % if bsxfun does not exist.
	bin = sum(repmat(smo_span,1,nspans-1) >= repmat(spans(1:end-1)',n,1),2);
end
dspans = diff(spans);
t = (smo_span - spans(bin))./dspans(bin);
index = (1:n)' + n*bin;
smo_raw = (1-t).*smo_n(index - n) + t.*smo_n(index);

% Apply final smooth.
smo = smooth(x,smo_raw,prop.weights,spans(1),prop.period);
smo = reshape(smo,sy);
if prop.unsorted
	smo(order) = smo;
end



%-------------------------------------------------------------------------
% Subfunctions
%-------------------------------------------------------------------------

function [smo,acvr] = smooth_wt_aper(x,y,w,span,period)
% smoothing function for aperiodic data, uses weights.

% Get the dimensions of the input.
n = length(y);

m = max(round(0.5*span*n),2);
k = 2*m + 1;

wy = w.*y;
wxy = wy.*x;
data = [cumprod([w,x,x],2),wy,wxy];

% Compute sum(w), sum(w*x), sum(w*x^2), sum(w*y), sum(w*x*y) over k points.
sums = zeros(n,5);
% ----------------------------------------------
% Slower, more accurate code:
% temp = filter(ones(it,1),1,data);
% sums(m+1:n-m,:) = temp(k:end,:);
% ----------------------------------------------
% Faster, slightly less accurate code:
cs = [0 0 0 0 0;cumsum(data)];
sums(m+1:n-m,:) = cs(k+1:end,:) - cs(1:end-k,:);
% ----------------------------------------------
sums(1:m,:) = sums((m+1)*ones(1,m),:);
sums(n-m+1:end,:) = sums((n-m)*ones(1,m),:);

denom = sums(:,1).*sums(:,3) - sums(:,2).^2;
a = (sums(:,4).*sums(:,3) - sums(:,2).*sums(:,5))./denom;
b = (sums(:,1).*sums(:,5) - sums(:,2).*sums(:,4))./denom;
smo = a + b.*x;

if nargout > 1
	sums_cv = sums - data;
	denom_cv = sums_cv(:,1).*sums_cv(:,3) - sums_cv(:,2).^2;
	a_cv = (sums_cv(:,4).*sums_cv(:,3) - sums_cv(:,2).*sums_cv(:,5))./denom_cv;
	b_cv = (sums_cv(:,1).*sums_cv(:,5) - sums_cv(:,2).*sums_cv(:,4))./denom_cv;
	smo_cv = a_cv + b_cv.*x;
	acvr = abs(smo_cv - y);
end

%-------------------------------------------------------------------------

function [smo,acvr] = smooth_wt_per(x,y,w,span,period)
% smoothing function for periodic data, uses weights.

% Get the dimensions of the input.
n = length(y);

m = max(round(0.5*span*n),2);
k = 2*m + 1;

x = [x;x(1:k-1)+period];
y = [y;y(1:k-1)];
w = [w;w(1:k-1)];

wy = w.*y;
wxy = wy.*x;
data = [cumprod([w,x,x],2),wy,wxy];

% Compute sum(w), sum(w*x), sum(w*x^2), sum(w*y), sum(w*x*y) over k points.
% ----------------------------------------------
% Slower, more accurate code:
% temp = filter(ones(k,1),1,data);
% sums = temp(k:end,:);
% ----------------------------------------------
% Faster, slightly less accurate code:
cs = [0 0 0 0 0;cumsum(data)];
sums = cs(k+1:n+k,:) - cs(1:n,:);
% ----------------------------------------------

denom = sums(:,1).*sums(:,3) - sums(:,2).^2;
a = (sums(:,4).*sums(:,3) - sums(:,2).*sums(:,5))./denom;
b = (sums(:,1).*sums(:,5) - sums(:,2).*sums(:,4))./denom;
% smo = circshift(a + b.*x(m+1:n+m),m); % slow
smo([m+1:n,1:m],:) = a + b.*x(m+1:n+m); % fast

if nargout > 1
	sums_cv = sums - data(m+1:n+m,:);
	denom_cv = sums_cv(:,1).*sums_cv(:,3) - sums_cv(:,2).^2;
	a_cv = (sums_cv(:,4).*sums_cv(:,3) - sums_cv(:,2).*sums_cv(:,5))./denom_cv;
	b_cv = (sums_cv(:,1).*sums_cv(:,5) - sums_cv(:,2).*sums_cv(:,4))./denom_cv;
% 	smo_cv = circshift(a_cv + b_cv.*x(m+1:n+m),m); % slow
	smo_cv([m+1:n,1:m],:) = a_cv + b_cv.*x(m+1:n+m); % fast
	acvr = abs(smo_cv - y(1:n));
end

%-------------------------------------------------------------------------

function [smo,acvr] = smooth_aper(x,y,w,span,period)
% smoothing function for aperiodic data, does not use weights.

% Get the dimensions of the input.
n = length(y);

m = max(round(0.5*span*n),2);
k = 2*m + 1;

xy = y.*x;
data = [cumprod([x,x],2),y,xy];

% Compute sum(x), sum(x^2), sum(y), sum(x*y) over k points.
sums = zeros(n,4);
% ----------------------------------------------
% Slower, more accurate code:
% temp = filter(ones(k,1),1,data);
% sums(m+1:n-m,:) = temp(k:end,:);
% ----------------------------------------------
% Faster, slightly less accurate code:
cs = [0 0 0 0;cumsum(data)];
sums(m+1:n-m,:) = cs(k+1:end,:) - cs(1:end-k,:);
% ----------------------------------------------
sums(1:m,:) = sums((m+1)*ones(1,m),:);
sums(n-m+1:end,:) = sums((n-m)*ones(1,m),:);

denom = k.*sums(:,2) - sums(:,1).^2;
a = (sums(:,3).*sums(:,2) - sums(:,1).*sums(:,4))./denom;
b = (k.*sums(:,4) - sums(:,1).*sums(:,3))./denom;
smo = a + b.*x;

if nargout > 1
	sums_cv = sums - data;
	denom_cv = (k-1).*sums_cv(:,2) - sums_cv(:,1).^2;
	a_cv = (sums_cv(:,3).*sums_cv(:,2) - sums_cv(:,1).*sums_cv(:,4))./denom_cv;
	b_cv = ((k-1).*sums_cv(:,4) - sums_cv(:,1).*sums_cv(:,3))./denom_cv;
	smo_cv = a_cv + b_cv.*x;
	acvr = abs(smo_cv - y);
end

%-------------------------------------------------------------------------

function [smo,acvr] = smooth_per(x,y,w,span,period)
% smoothing function for periodic data, does not use weights.

% Get the dimensions of the input.
n = length(y);

m = max(round(0.5*span*n),2);
k = 2*m + 1;

x = [x;x(1:k-1)+period];
y = [y;y(1:k-1)];

xy = y.*x;
data = [cumprod([x,x],2),y,xy];

% Compute sum(x), sum(x^2), sum(y), sum(x*y) over k points.
% ----------------------------------------------
% Slower, more accurate code:
% temp = filter(ones(k,1),1,data);
% sums = temp(k:end,:);
% ----------------------------------------------
% Faster, slightly less accurate code:
cs = [0 0 0 0;cumsum(data)];
sums = cs(k+1:n+k,:) - cs(1:n,:);
% ----------------------------------------------

denom = k.*sums(:,2) - sums(:,1).^2;
a = (sums(:,3).*sums(:,2) - sums(:,1).*sums(:,4))./denom;
b = (k.*sums(:,4) - sums(:,1).*sums(:,3))./denom;
% smo = circshift(a + b.*x(m+1:n+m),m); % slow
smo([m+1:n,1:m],:) = a + b.*x(m+1:n+m); % fast

if nargout > 1
	sums_cv = sums - data(m+1:n+m,:);
	denom_cv = (k-1).*sums_cv(:,2) - sums_cv(:,1).^2;
	a_cv = (sums_cv(:,3).*sums_cv(:,2) - sums_cv(:,1).*sums_cv(:,4))./denom_cv;
	b_cv = ((k-1).*sums_cv(:,4) - sums_cv(:,1).*sums_cv(:,3))./denom_cv;
% 	smo_cv = circshift(a_cv + b_cv.*x(m+1:n+m),m); % slow
	smo_cv([m+1:n,1:m],:) = a_cv + b_cv.*x(m+1:n+m); % fast
	acvr = abs(smo_cv - y(1:n));
end
