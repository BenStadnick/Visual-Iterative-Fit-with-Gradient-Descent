function [  ] = VisualIterativeFit(  )
%VISUALITERRATIVEFIT iterrative fit algorithm using the gradient descnent
% method that plots the progress as the fit occurs to illastrate the
% numerical fittng process and effects of error on the data. Data sets are
% generated using the fitting function with added error to the data.
%   Created by Ben Stadnick on June 2017

    % Set up variables for generating data
    DataDetLenght = 100;
    trueparameters = [-7; 1];% True parameters for fitting function
    RandomError = 3;
    SystematicError = 1;
    DataRange = [0; 10];
    
    % Generate a data set to be fitted, 'x' is our independent variable and
    % 'y' is our ddependent variable
    x = linspace(DataRange(1), DataRange(2), DataDetLenght)';
    y = fitfunction(x, trueparameters);
    % Apply error
    for n = 1:DataDetLenght
        y(n) = y(n) + RandomError*rand(1) - RandomError/2 + SystematicError;
    end

    % Set up variables for fitting algorithm
    FitParameters = [-1; 0.5];% Initial guesses for fitted parameters
    StepSize = [1, 1];% Step size for each parameter
    Iterations = 200;% Maximim number of iterations
    Tolerance = 100;% Acceptable sum of least square errors
    fig = figure(1);
    
    % Run itterate fit and plot
    for n = 1:Iterations
        % Calculate the sum of least squares
        LeastSquaresValue_Current = sum( (fitfunction(x, FitParameters) - y).^2 );

        parameters_new = FitParameters;
        % Running gradient descent
        for k = 1:length(FitParameters)
            % Calculate new value
            parameters_new(k) = FitParameters(k)+StepSize(k);
            y_trial = fitfunction(x, FitParameters);
            LeastSquaresValue_New = sum( (fitfunction(x, parameters_new) - y).^2 );

            % Updating the parameters
            if(LeastSquaresValue_New < LeastSquaresValue_Current)
                LeastSquaresValue_Current = LeastSquaresValue_New;
                FitParameters = parameters_new;
                continue
            end

            % Calculate new value
            parameters_new(k) = FitParameters(k)-StepSize(k);
            y_trial = fitfunction(x, FitParameters);
            LeastSquaresValue_New = sum( (fitfunction(x, parameters_new) - y).^2 );

            % Updating the parameters
            if(LeastSquaresValue_New < LeastSquaresValue_Current)
                LeastSquaresValue_Current = LeastSquaresValue_New;
                StepSize(k) = -StepSize(k);
                FitParameters = parameters_new;
                continue
            end

            % If niether directions reduce the sum of least squares,
            % decrease step size
            StepSize(k) = StepSize(k)/2; 
        end
        
        % Plotting our final hypothesis
        figure(1);
        scatter(x, y);
        hold on;
        plot(x, fitfunction( x, FitParameters ), 'linewidth', 2, 'color', 'red');
        hold off
        LeastSquaresValue = sum( (fitfunction(x, FitParameters) - y).^2);
        title( sprintf('Itteration: %.0f    SumOfSquaresError: %.2f', n, LeastSquaresValue));

        % If tolerance is reached, end fit
        if(LeastSquaresValue_New <= Tolerance)
            break;
        end
    end
end
%%
%Fitting function
function [ y ] = fitfunction( x, parameters )
    y = x*parameters(1) + x.^2*parameters(2);    
end
