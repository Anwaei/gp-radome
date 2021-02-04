function hyp=function_gp_hyperparam(Xm,fmh)

nm = size(Xm,2); % This is the number of measurement points.

%% Next, we set hyperparameters to blatantly incorrect values, so we can tune them.
lf=1;
lx=1;
sfh=0.1;
hyp = [lx^2;lf^2;sfh^2]; % This is an array of hyperparameters which will be tuned.

%% We set things up for the gradient ascent algorithm.
numSteps = 100; % How many gradient ascent steps shall we take?
stepSize = 1; % This is the initial step size. In the non-dimensionalized gradient ascent algorithm we use below, this can be seen as a length scale of the optimized parameter, in this case the log-likelihood.
stepSizeFactor = 2; % This is the factor by which we will decrease the step size in case it is too big.
maxReductions = 100; % This is the maximum number of times in a row which we can reduce the step size. If we'd need this many reductions, something is obviously wrong.
clear logp; % We make sure that logp is not defined. Whether it is defined is used in the script to check if it's the first run of the algorithm.
newHypDeriv = zeros(3,1); % We already create a storage for the new hyperparameter derivative array. We'll need this soon.

%% Now we can start iterating
for i = 1:numSteps
	% We try to improve the parameters, all the while checking the step size.
	for j = 1:maxReductions
		% We check if we haven't accidentally been decreasing the step size too much.
		if j == maxReductions
			disp('Error: something is wrong with the step size in the hyperparameter optimization scheme.');
		end
		% We calculate new hyperparameters. Or at least, candidates. We still check them.
		if ~exist('logp','var') % If no logp is defined, this is the first time we are looping. In this case, with no derivative data known yet either, we keep the hyperparameters the same.
			newHyp = hyp;
		else
			newHyp = hyp.*(1 + stepSize*hyp.*hypDeriv); % We apply a non-dimensional update of the hyperparameters. This only works when the parameters are always positive.
		end
		% Now we check the new hyperparameters. If they are good, we will implement them.
		if min(newHyp > 0) % The parameters have to remain positive. If they are not, something is wrong. To be precise, the step size is too big.
			% We partly implement the new hyperparameters and check the new value of logp.
			lx = sqrt(newHyp(1));
			lf = sqrt(newHyp(2));
			sfh = sqrt(newHyp(3));
			Sfm = sfh^2*eye(nm); % This is the noise covariance matrix.
			diff = repmat(Xm,nm,1) - repmat(Xm',1,nm); % This is a matrix with element [i,j] equal to x_j - x_i.
			Kmm = lf^2*exp(-1/2*diff.^2/lx^2); % This is the covariance matrix. It contains the covariances of each combination of points.
			P = Kmm + Sfm;
 			mb = (ones(nm,1)'/P*fmh)/(ones(nm,1)'/P*ones(nm,1)); % This is the (constant) mean function m(x) = \bar{m}. You can get rid of this line if you don't want to tune \bar{m}.
			newLogp = -nm/2*log(2*pi) - 1/2*logdet(P) - 1/2*(fmh - mb)'/P*(fmh - mb);
			% If this is the first time we are in this loop, or if the new logp is better than the old one, we fully implement the new hyperparameters and recalculate the derivative.
			if ~exist('logp','var') || newLogp >= logp
				% We calculate the new hyperparameter derivative.
				alpha = P\(fmh - mb);
				R = alpha*alpha' - inv(P);
				newHypDeriv(3) = 1/2*trace(R);
				newHypDeriv(2) = 1/(2*lf^2)*trace(R*Kmm);
				newHypDeriv(1) = 1/(4*lx^4)*trace(R*(Kmm.*(diff.^2)));
				% If this is not the first time we run this, we also update the step size, based on how much the (normalized) derivative direction has changed. If the derivative is still in the
				% same direction as earlier, we take a bigger step size. If the derivative is in the opposite direction, we take a smaller step size. And if the derivative is perpendicular to
				% what is used to be, then the step size was perfect and we keep it. For this scheme, we use the dot product.
				if exist('logp','var')
					directionConsistency = ((hypDeriv.*newHyp)'*(newHypDeriv.*newHyp))/norm(hypDeriv.*newHyp)/norm(newHypDeriv.*newHyp);
					stepSize = stepSize*stepSizeFactor^directionConsistency;
				end
				break; % We exit the step-size-reduction loop.
			end
		end
		% If we reach this, it means the hyperparameters we tried were not suitable. In this case, we should reduce the step size and try again. If the step size is small enough, there will
		% always be an improvement of the hyperparameters. (Unless they are fully perfect, which never really occurs.)
		stepSize = stepSize/stepSizeFactor;
	end
	% We update the important parameters.
	hyp = newHyp;
	hypDeriv = newHypDeriv;
	logp = newLogp;

end