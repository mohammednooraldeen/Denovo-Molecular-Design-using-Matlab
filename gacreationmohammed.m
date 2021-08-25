function Population = gacreationmohammed(GenomeLength,FitnessFcn,options)
%GACREATIONUNIFORM Creates the initial population for genetic algorithm.
%   POP = GACREATIONUNIFORM(NVARS,FITNESSFCN,OPTIONS) Creates the
%   initial population that GA will then evolve into a solution.
%
%   Example:
%     options = gaoptimset('CreationFcn',@gacreationuniform);

%   Copyright 2003-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2010/11/01 19:36:46 $
%

if strcmpi(options.PopulationType,'custom')
    error(message('globaloptim:gacreationuniform:unknownPopulationType', options.PopulationType));
end

totalPopulation = sum(options.PopulationSize);
initPopProvided = size(options.InitialPopulation,1);
individualsToCreate = totalPopulation - initPopProvided;

if strcmpi(options.PopulationType,'doubleVector')
    % Initialize Population to be created
    Population = zeros(totalPopulation,GenomeLength);
    % Use initial population provided already
    if initPopProvided > 0
        Population(1:initPopProvided,:) = options.InitialPopulation;
    end
    % Create remaining population
    % problemtype is either 'unconstrained', 'boundconstraints', or
    % 'linearconstraints'. Nonlinear constrained algorithm 'ALGA' does not
    % create or use initial population of its own. It calls sub-problem
    % solvers (galincon/gaunc)
    if isfield(options,'LinearConstr')
        problemtype = options.LinearConstr.type;
    else
        problemtype = 'unconstrained';
    end
    % This function knows how to create initial population for
    % unconstrained and boundconstrained cases but does not know in
    % linearconstrained case.
    if ~strcmp(problemtype,'linearconstraints')
        range = options.PopInitRange;
        lowerBound = range(1,:);
        span = range(2,:) - lowerBound;
        %Population(initPopProvided+1:end,:) = repmat(lowerBound,individualsToCreate,1) + ...
            repmat(span,individualsToCreate,1) .* rand(individualsToCreate,GenomeLength);
     p=zeros[1,GenomeLength];
     % setting first point
     % each structure is being made which will be part of the full gene
     % the atom vector includes [no. x y z no.ofboundsavailable
     % whichcanbeused coordinatesofavailabebonds]
     
     % the function to calculate rest the two connection points of an atom
     % is:
     coord_of_first
     coord_of_second
     x_deviation=cosd(30)*(2*sind(60)*1.33);
     y_deviation=sind(30)*(2*sind(60)*1.33);
     %rotation matrix about arbitrary axis
     c=cosd(theta);s=sind(theta);
     
     
     
     for i=1:num_atoms
     if (i==1) p(1,1:5)=[1 0 0 0 3];end
   
     n_ready=p(:,5)>=1;
     s=size(n_ready);
     choice=randsample(n_ready(:,1),1);
     p(choice,5)=p(choice,5)-1;
     pole=p(choice,2:4);
     %generate random angle sign
     theta=120;
     if(rand(1)>0.5) theta_sign=-1 ;else theta_sign=1; end
     
     point=rotmat([0;0;1],theta_sign*theta,[1.33,0,0]'); point=point';
     
    for j=1:s(1,1); if (p(j,2:4)==point) j=s(1,1); else p(i,1:5)=[i,point,2];end;end
    
     
     
    
    
    else
        
        Population(initPopProvided+1:end,:) = []; 
    end
elseif strcmpi(options.PopulationType,'bitString')
    % Initialize Population to be created
    Population = ones(totalPopulation,GenomeLength);
    % Use initial population provided already
    if initPopProvided > 0
        Population(1:initPopProvided,:) = options.InitialPopulation;
    end
    % Create remaining population
    Population(initPopProvided+1:end,:) = double(0.5 > rand(individualsToCreate,GenomeLength));
end

if all(isnan(Population))
    error(message('globaloptim:gacreationuniform:populationIsNaN'));
elseif all(isinf(Population))
    error(message('globaloptim:gacreationuniform:populationIsInf'));
end

