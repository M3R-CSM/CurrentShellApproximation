function [ currentEstimate, cost, soln_found ] = LM_LeastSquares( paramIC, costFunction, tol, maxIterations, maxTime, dispVals, maxChange )
%LM_LEASTSQUARES solves a minimization problem using the LM algorithm.
%   [ currentEstimate, cost ] = LM_LeastSquares( paramIC, costFunction, tol, dispVals )
%
%   currentEstimate: the final parameter estimate
%   cost:  the  L2 Norm of the residual error
%
%   paramIC: An initial guess at the solution
%   costFunction: a function that returns an error vector and a jacobian
%       [error, Jacobian] = function(parameters)
%   tol: the exit tolerance default is 1e-4
%   maxIterations: the maximum number of iterations
%   maxTime: the maximum amount of time to run
%   dispVals:  true to print progress, false to be silent
%
% Author: Andrew J Petruska
%

TSTART=tic;
soln_found = true;
if ~exist('dispVals','var')
    dispVals = false;
end

if ~exist('tol','var')
    tol = 1e-4;
end

if ~exist('maxIterations','var')
    maxIterations = 10000;
end

if ~exist('maxTime','var')
    maxTime = inf;
end

if ~exist('maxChange','var')
    maxChange = inf;
end


currentEstimate = paramIC;
last_estimateChange = zeros(size(paramIC));
[errorVector,J] = costFunction(currentEstimate);
errorVectorLast = errorVector;


converged = false;

errorThis = norm(errorVector);
errorLast = errorThis;
errorPrediction = -1;

if dispVals
    optimValues.iteration =0;
    optimValues.funccount = 0;
    optimValues.fval = errorThis;
    optimValues.procedure = '';
    callAllOptimPlotFcns({@optimplotfval},currentEstimate,optimValues,'init');
end

Lambda = 1;

LamdaChangeRange = .5;

count = 0;
tmp = warning('error','MATLAB:lscov:RankDefDesignMat');
proc = [];
lastBadStep = -1;
while (~converged && count < maxIterations && (toc(TSTART) < maxTime) )
    
    deltaError = errorLast-errorThis;
    predictionFactor = deltaError/(errorLast-errorPrediction);
    if errorPrediction == -1 || isinf(predictionFactor) || isnan(predictionFactor)
        % no change in error or initial pass
        predictionFactor = 1-3/4*LamdaChangeRange;
    end
    
    
    if deltaError < 0 %% Error Got worse
        Lambda = 2*Lambda;
        if count - lastBadStep > 1
            LamdaChangeRange = max(LamdaChangeRange/2,tol/10);
            if LamdaChangeRange == tol/100
                warning('LM_LeastSquares: No Solution Found');
                soln_found = false;
                break;
            end
        else
            S = svd(J);
            S1 = S(1);
            Se = min(S(S~=0));
            if Lambda*(500*Se/S1) < 1
                Lambda = min(10*Lambda,1/(500*Se/S1));
            end
        end
        lastBadStep = count;
        
        estimateChange = -last_estimateChange;
        last_estimateChange = zeros(size(last_estimateChange));
        proc = 'Increase Lambda';
    else % Error improved
        
        %Calculate new weighted jacobian
        H = zeros(size(J,2));
        Hinv = zeros(size(J,2));
        for i=1:size(J,2)
            Hinv(i,i) = norm(J(:,i));
            if Hinv(i,i) > 10*eps
                H(i,i) = 1/Hinv(i,i);
                J(:,i) = J(:,i)*H(i,i);
            else
                H(i,i) = 0;
                Hinv(i,i) = 0;
                J(:,i) = J(:,i)*H(i,i);
                disp('Warning Matrix is close to singular');
            end
        end
        
        % See how linearly the system behavied (prediction factor =1) is a
        % linear system
        if abs(1-predictionFactor) > LamdaChangeRange && abs(deltaError) > sqrt(eps)
            % Error change was too small compared to a linear change.
            % Getting close to an over-step.  Be more gradient decent.
            Lambda = Lambda*2;
            
            %estimateChange = -H*lscov(J'*J+Lambda*diag(diag(J'*J)),J'*errorVector);
            try
                estimateChange = -H*lscov(J'*J+Lambda*diag(diag(J'*J)),J'*errorVector);
                if norm(estimateChange) > maxChange 
                    estimateChange = estimateChange*maxChange/norm(estimateChange);
                end
            catch exception
                Lambda = Lambda*100;
            end
            last_estimateChange = estimateChange;
            
            proc = 'Step But Increase Lambda';
            
        elseif isempty(deltaError) || ( abs(1-predictionFactor) < LamdaChangeRange/2)||(abs(deltaError) >0 && abs(deltaError)<sqrt(eps))% && abs(deltaError) > sqrt(eps))%state == 2
            % Error change was too close to linear.  Get more aggresive!
            Lambda = Lambda/2;
            try
                estimateChange = -H*lscov(J'*J+Lambda*diag(diag(J'*J)),J'*errorVector);
                if norm(estimateChange) > maxChange
                    estimateChange = estimateChange*maxChange/norm(estimateChange);
                end
                last_estimateChange = estimateChange;
                
            catch exception
                disp('Warning Matrix is close to singular');
                estimateChange = -H*pseudoInv(J'*J+Lambda*diag(diag(J'*J)),1/100000)*J'*errorVector;
                if norm(estimateChange) > maxChange
                    estimateChange = estimateChange*maxChange/norm(estimateChange);
                end
                last_estimateChange = estimateChange;
            end
            
            proc = 'Step and Decrease Lambda';
            
        else
            % Error change was in a sweet spot.  Keep on keeping on.
            try
                estimateChange = -H*lscov(J'*J+Lambda*diag(diag(J'*J)),J'*errorVector);
                if norm(estimateChange) > maxChange 
                    estimateChange = estimateChange*maxChange/norm(estimateChange);
                end
            catch exception
                Lambda = Lambda*100;
            end
            last_estimateChange = estimateChange;
            if abs(deltaError) > sqrt(eps)
                LamdaChangeRange = LamdaChangeRange*2;
            end
            proc = 'Step Maintain Lambda';
            
            
        end
    end
    
    %% update Estimate
    currentEstimate = currentEstimate + estimateChange;
    
    if deltaError >= 0
        %update errorVector
        errorPrediction = norm(errorVector+J*(Hinv*estimateChange));
        errorLast = errorThis;
        errorVectorLast = errorVector;
        [errorVector,J] = costFunction(currentEstimate);
        if any(isnan(errorVector)) || any(any(isnan(J)))
            currentEstimate = currentEstimate-estimateChange;
            warning('LM_LeastSquares:costFunction returned NaN for the error or the jacobian');
            break;
        end
        errorThis = norm(errorVector);
        
        converged = errorThis < errorLast || errorThis < tol;
        %         if converged
        %             converged = converged && norm(estimateChange)/norm(currentEstimate)*(1+Lambda) < tol;
        %         end
        if converged
            converged = (abs(errorThis-errorLast)/errorThis*(1+Lambda) < tol && norm(estimateChange)/norm(currentEstimate) < tol || errorThis < tol) && count > 0;
        end
        count = count +1;
    else
        % Reset to last values
        errorVector = errorVectorLast;
        errorThis = errorLast;
        errorPrediction = -1;
        
        converged = false;
        
        if Lambda > 1/(10*tol)
            converged = true; %% give up
        end
    end
    
    
    if dispVals
        optimValues.iteration = count-1;
        optimValues.funccount = count;
        optimValues.fval = errorThis;
        optimValues.procedure = proc;
        
        stop = callAllOptimPlotFcns({@optimplotfval},currentEstimate,optimValues,'iter');
        if stop
            warning('LM_LeastSquares: Convergance Terminated by user');
                soln_found = false;

            break;
        end
        converged = converged && ~stop;
        fprintf(['Iteration: %i\tCost: %5.5f\tConverged: %i\tError Change: %5.5f%%\tValue Change: %5.5f%%\tLambda: %5.5f\tLambda Change Range: %5.5d\t' proc '\n'],[count,errorThis,converged,(errorThis-errorLast)/errorThis*100,norm(estimateChange)/norm(currentEstimate)*100,Lambda,LamdaChangeRange]);
    end
    
end

if toc(TSTART) > maxTime
    warning('LM_LeastSquares:Maximum Time Exceeded');
    soln_found = false;
end

if count >= maxIterations
    warning('LM_LeastSquares:Maximum Iterations Exceeded');
        soln_found = false;

end

if (abs(errorThis-errorLast)/errorThis*(1+Lambda) < tol && norm(estimateChange)/norm(currentEstimate) < tol)
    warning('LM_LeastSquares:Error or Esitimate not changing. Local minimum found.');
        soln_found = false;

end

if dispVals
    callAllOptimPlotFcns({@optimplotfval},currentEstimate,optimValues,'done');
end
cost = errorThis;

end

