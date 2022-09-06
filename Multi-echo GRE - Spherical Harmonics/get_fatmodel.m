function fatmodel = get_fatmodel(num_peaks)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% get_fatmodel: function that selects the correct fat model, where a struct 
% is returned that contains both the frequencies of the fat peaks (with 
% respect to the water frequency) as well as their relative amplitudes 
% (which sum to 1)
% The reference from which the models are determined are provided for each
% model separately
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   num_peaks:   integer value of the number of fat peaks
%
% Outputs:
%   fatmodel:    struct containing frequencies and relative amplitudes of
%                the fat peaks of the desired model


if num_peaks == 3
    % 3 peak model
    % Yu H: Multiecho water-fat separation and simultaneous R2* estimation with multifrequency fat spectrum modeling
    fatmodel.frequency = [0.73 -2.49 -3.29];
    fatmodel.relAmps = [0.08 0.17 0.75];
    
elseif num_peaks == 4
    % 4 peak model
    % Wokke BH: Comparison of dixon and T1-weighted MR methods to assess the degree of fat infiltration in duchenne muscular dystrophy patients
    fatmodel.frequency = [0.73 -2.49 -3.27 -3.68];
    fatmodel.relAmps = [0.08 0.15 0.72 0.04];
    
elseif num_peaks == 5
    % 5 peak model
    % Wokke BH: Comparison of dixon and T1-weighted MR methods to assess the degree of fat infiltration in duchenne muscular dystrophy patients
    fatmodel.frequency = [0.73 -2.35 -2.54 -3.27 -3.68];
    fatmodel.relAmps = [0.08 0.05 0.1 0.72 0.4];
    
elseif num_peaks == 6
    % 6 peak model
    % Hernando D: Chemical shift-based water/fat separation: a comparison of signal models
    fatmodel.frequency = [0.6 -0.5 -1.95 -2.6 -3.4 -3.8];
    fatmodel.relAmps = [0.047, 0.039, 0.006, 0.12, 0.7, 0.088];

elseif num_peaks == 7
    % 7 peak model
    % Ren J: Composition of adipose tissue and marrow fat in humans by 1H NMR at 7 Tesla
    fatmodel.frequency = [0.61 -1.93 -2.45 -2.67 -3.11 -3.4 -3.8];
    fatmodel.relAmps = [0.042 0.015 0.066 0.096 0.071 0.627 0.083];
    
elseif num_peaks == 9
    % 9 peak model
    % Hamilton G: In vivo characterization of the liver fat (1)H MR spectrum
    fatmodel.frequency = [0.59 0.49 -0.5 -1.95 -2.46 -2.68 -3.10 -3.40 -3.80];
    fatmodel.relAmps = [0.0538 0.01 0.0399 0.014 0.0598 0.0798 0.0598 0.5932 0.0897];
    
else
    message = {'Selection of the correct fat model was not possible, '
        'please set num_fat_peaks to the correct integer value. '};
    errorStruct.message = sprintf('%s\n',message{:});
    errorStruct.identifier = 'MyFunction:parameterIncorrecctlyDefined';
    error(errorStruct)
    
end
end