necessaryFeatures = {'Optimization_Toolbox','GADS_Toolbox','Statistics_Toolbox'};
index = cellfun(@(f) license('checkout',f),necessaryFeatures);
% Copyright 2014-2015: Lars Freier, Eric von Lieres
% See the license note at the end of the file.

for iToolbox=1:length(index)
    if index(iToolbox)<1
        warning('No active licence for %s was found',necessaryFeatures{iToolbox}) 
        disp('Press any key')
        pause
    end
end

disp('Checked Licenses')

%%
% featureStr = {'Aerospace_Blockset'; ...
%               'Aerospace_Toolbox'; ...
%               'Bioinformatics_Toolbox'; ...
%               'Communication_Blocks'; ...
%               'Communication_Toolbox'; ...
%               'Compiler'; ...
%               'Control_Toolbox'; ...
%               'Curve_Fitting_Toolbox'; ...
%               'Data_Acq_Toolbox'; ...
%               'Database_Toolbox'; ...
%               'Datafeed_Toolbox'; ...
%               'Dial_and_Gauge_Blocks'; ...
%               'Distrib_Computing_Toolbox'; ...
%               'Econometrics_Toolbox'; ...
%               'EDA_Simulator_Link_DS'; ...
%               'Embedded_Target_c166'; ...
%               'Embedded_Target_c2000'; ...
%               'Embedded_Target_c6000'; ...
%               'Embedded_Target_MPC555'; ...
%               'Excel_Link'; ...
%               'Filter_Design_HDL_Coder'; ...
%               'Filter_Design_Toolbox'; ...
%               'Fin_Derivatives_Toolbox'; ...
%               'Financial_Toolbox'; ...
%               'Fixed_Income_Toolbox'; ...
%               'Fixed_Point_Toolbox'; ...
%               'Fixed-Point_Blocks'; ...
%               'Fuzzy_Toolbox'; ...
%               'GADS_Toolbox'; ...
%               'IDE_Link_MU'; ...
%               'Identification_Toolbox'; ...
%               'Image_Acquisition_Toolbox'; ...
%               'Image_Toolbox'; ...
%               'Instr_Control_Toolbox'; ...
%               'Link_for_Incisive'; ...
%               'Link_for_ModelSim'; ...
%               'Link_for_Tasking'; ...
%               'Link_for_VisualDSP'; ...
%               'MAP_Toolbox'; ...
%               'MATLAB'; ...
%               'MATLAB_Builder_for_dot_Net'; ...
%               'MATLAB_Builder_for_Java'; ...
%               'MATLAB_Distrib_Comp_Engine'; ...
%               'MATLAB_Excel_Builder'; ...
%               'MATLAB_Link_for_CCS'; ...
%               'MATLAB_Report_Gen'; ...
%               'MBC_Toolbox'; ...
%               'MPC_Toolbox'; ...
%               'NCD_Toolbox'; ...
%               'Neural_Network_Toolbox'; ...
%               'OPC_Toolbox'; ...
%               'Optimization_Toolbox'; ...
%               'PDE_Toolbox'; ...
%               'Power_System_Blocks'; ...
%               'Real-Time_Win_Target'; ...
%               'Real-Time_Workshop'; ...
%               'RF_Blockset'; ...
%               'RF_Toolbox'; ...
%               'Robust_Toolbox'; ...
%               'RTW_Embedded_Coder'; ...
%               'Signal_Blocks'; ...
%               'Signal_Toolbox'; ...
%               'SimBiology'; ...
%               'SimDriveline'; ...
%               'SimElectronics'; ...
%               'SimEvents'; ...
%               'SimHydraulics'; ...
%               'SimMechanics'; ...
%               'Simscape'; ...
%               'SIMULINK'; ...
%               'Simulink_Control_Design'; ...
%               'Simulink_Design_Verifier'; ...
%               'Simulink_HDL_Coder'; ...
%               'Simulink_Param_Estimation'; ...
%               'SIMULINK_Report_Gen'; ...
%               'SL_Verification_Validation'; ...
%               'Spline_Toolbox'; ...
%               'Stateflow'; ...
%               'Stateflow_Coder'; ...
%               'Statistics_Toolbox'; ...
%               'Symbolic_Toolbox'; ...
%               'SystemTest'; ...
%               'Video_and_Image_Blockset'; ...
%               'Virtual_Reality_Toolbox'; ...
%               'Wavelet_Toolbox'; ...
%               'XPC_Embedded_Option'; ...
%               'XPC_Target'};
% index = cellfun(@(f) license('checkout',f),featureStr);
% featureStr(logical(index));
% =============================================================================
%  KriKit - Kriging toolKit
%  
%  Copyright 2014-2015: Lars Freier(1), Eric von Lieres(1)
%                                      
%    (1)Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
