clc;
clear;

timestep = 0.001;
num_steps = 5000;
cell_size = 2;
dim = 3; %Spin-one model
chi = 100;

state = SpinMPS.initialize_ferromagnet(cell_size,dim);
[operator,hamilt] = MPS_BondOperator.MPS_AKLT_Bond(cell_size);

for j = 1:num_steps
    state = operator.apply_to_state(state,chi,timestep);
    energies = state.bond_expectation(hamilt);
    mean_E = mean(cell2mat(energies));
    fprintf('Iteration %i: energy = %d \n',j,mean_E);
end

%% Spectrum

ent_spec = state.entanglement_spectrum(8)

