if exist('figure_handles','var') 
    for j = 1:numel(figure_handles)
        if ishandle(figure_handles{j})
            close(figure_handles{j});
        end
    end
    clear('figure_handles');
end

clc;
clear;

figure_handles = cell(1,1);

addpath(fullfile(pwd,'..'));

%********INPUTS**********
field = 0.5;
coupling = 1;
cell_size = 2;
timestep = 0.025;
num_steps = 200;
chi = 50;
%************************

dim = 2; %Spin-half model

state = SpinMPS.initialize_ferromagnet(cell_size,dim);
[operator,hamilt] = MPS_BondOperator.MPS_Ising_Bond(coupling,field,cell_size);

for j = 1:num_steps
    state = operator.apply_to_state(state,chi,timestep);
    energies = state.bond_expectation(hamilt);
    mean_E = mean(cell2mat(energies));
    fprintf('Iteration %i: energy = %d \n',j,mean_E);
end

sing_en = @(k) -2*sqrt(1 + field^2 - 2*field*cos(k));
exact_en = integral(sing_en,0,pi)/(2*pi);

fprintf('Exact energy is %d',exact_en);
