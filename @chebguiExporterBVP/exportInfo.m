function expInfo = exportInfo(e, guifile)

% Extract information from the GUI fields
dom = guifile.domain;
deInput = guifile.DE;
bcInput = guifile.BC;
initInput = guifile.init;

% Wrap all input strings in a cell (if they're not a cell already)
if ( isa(deInput, 'char') )
     deInput = cellstr(deInput);
end

if ( isa(bcInput, 'char') )
     bcInput = cellstr(bcInput);
end

if ( isa(initInput, 'char') )
     initInput = cellstr(initInput);
end

% Obtain useful strings describing the differential equation part:
[deString, allVarString, indVarNameDE, dummy, dummy, dummy, allVarNames] = ...
    setupFields(guifile, deInput, 'DE');

% Do some error checking before we do further printing. Check that
% independent variable name match.
% Obtain the independent variable name appearing in the initial condition
useLatest = strcmpi(initInput{1}, 'Using latest solution');
if ( ~isempty(initInput{1}) && ~useLatest )
    [dummy, dummy, indVarNameInit] = ...
        setupFields(guifile, initInput, 'BC', allVarString);
else
    indVarNameInit = {''};
end

% Make sure we don't have a discrepency in indVarNames
if ( ~isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1}) )
    if ( strcmp(indVarNameDE{1}, indVarNameInit{1}) )
        indVarNameSpace = indVarNameDE{1};
    else
        error('Chebgui:SolveGUIbvp', 'Independent variable names do not agree')
    end
elseif ( ~isempty(indVarNameInit{1}) && isempty(indVarNameDE{1}) )
    indVarNameSpace = indVarNameInit{1};
elseif ( isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1}) )
    indVarNameSpace = indVarNameDE{1};
else
    indVarNameSpace = 'x'; % Default value
end

% Replace the 'DUMMYSPACE' variable in the DE field
deString = strrep(deString, 'DUMMYSPACE', indVarNameSpace);
deString = chebguiExporter.prettyPrintFevalString(deString, allVarNames);

% Support for periodic boundary conditions
if ( ~isempty(bcInput{1}) && strcmpi(bcInput{1}, 'periodic') )
    bcInput{1} = [];
    periodic = true;
else
    periodic = false;
end

%% Information related to options set-up:
% Option for tolerance
tolInput = guifile.tol;


%% Fill up the expInfo struct
expInfo.dom = dom;
expInfo.deInput = deInput;
expInfo.bcInput = bcInput;
expInfo.initInput = initInput;
expInfo.deString = deString;
expInfo.allVarString = allVarString;
expInfo.allVarNames = allVarNames;
expInfo.indVarNameSpace = indVarNameSpace;
expInfo.periodic = periodic;
expInfo.useLatest = useLatest;

% Information related to options set-up
expInfo.tol = guifile.tol;
expInfo.dampedOn = guifile.options.damping;
expInfo.discretization = guifile.options.discretization;
expInfo.plotting = guifile.options.plotting;
end