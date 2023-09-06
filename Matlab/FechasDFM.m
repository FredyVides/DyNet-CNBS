%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\fredy.vides\OneDrive - Comisión Nacional de Bancos y Seguros (CNBS)\DatosEF\DatosDinamicaEconomica.xlsx
%    Worksheet: Resumen Variables
%
% Auto-generated by MATLAB on 28-Jul-2023 15:00:56

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 1);

% Specify sheet and range
opts.Sheet = "Resumen Variables";
opts.DataRange = "A5:A245";

% Specify column names and types
opts.VariableNames = "VarName1";
opts.VariableTypes = "datetime";

% Import the data
DatesDFM = readtable("DatosDinamicaEconomica.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts