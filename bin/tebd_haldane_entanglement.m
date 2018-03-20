clc;
clear;

timestep = 0.05;
num_steps = 1000;
cell_size = 2;
dim = 3; %Spin-one model
chi = 80;

converge_limit = 1.e-5; converge_range = 50;

%********Parameters*********
zz_1 = 0.1;
zz_2 = -0.1;
heis = 1;
inversion = 0.2;
rotation = 0.2;
xfield_1 = 0.05;
xfield_2 = 0.0;
zfield = 0.05;
time_lim = 10;
%*************************

RUN = true;

if RUN

    state = SpinMPS.initialize_ferromagnet(cell_size,dim);
    [operator_1,hamilt_1] = MPS_BondOperator.MPS_Haldane_Bond(heis,zz_1,inversion,rotation,xfield_1,zfield,cell_size);

    iter_energies = zeros(1,num_steps);

    for j = 1:num_steps
        [state,max_error] = operator_1.apply_to_state(state,chi,timestep);
        energies = state.bond_expectation(hamilt_1);
        iter_energies(j) = mean(cell2mat(energies));
        evals = real(state.eigs{cell_size});
        [e_gap,ind] = max(evals); evals(ind) = []; e_gap = e_gap - max(evals);
        fprintf('Iteration %i: energy = %d; chi = %i; max_error = %d; ent_gap = %d \n',j,iter_energies(j),numel(state.eigs{1}),max_error,e_gap);
        if j > (converge_range+1)
            if abs(iter_energies(j) - iter_energies(j-converge_range)) < converge_limit
                fprintf('\n**************\nCONVERGED at iteration %i: energy = %d \n\n',j,iter_energies(j));
                break;
            end
        end
    end

    save('mps_state.mat','state','heis','zz_1','inversion','rotation','xfield_1');
    
else
    load('mps_state.mat','state','heis','zz_1','inversion','rotation','xfield_1');
end

big_es = state.eigs{2};
disp(big_es(1:10));

%% Time evolution

state.test_inversion_symmetry();

% real_time_steps = floor(time_lim/timestep);
% real_times = (0:real_time_steps).*timestep;
% 
% ent_spec = zeros(chi,1+real_time_steps);
% ent_spec(:,1) = -2*log(state.eigs{2});
% 
% [operator_2,hamilt_2] = MPS_BondOperator.MPS_Haldane_Bond(heis,zz_2,inversion,rotation,xfield_2,zfield,cell_size);
% 
% for j = 1:real_time_steps
%     [state,max_error] = operator_2.apply_to_state(state,chi,1i*timestep);
%     evals = real(state.eigs{cell_size});
%     ent_spec(:,j+1) = -2*log(sort(evals,'descend'));
%     e_gap = ent_spec(2,j+1) - ent_spec(1,j+1);
%     fprintf('Time %f: chi = %i; max_error = %d; ent_gap = %d \n',real_times(j),numel(state.eigs{1}),max_error,e_gap);
% end
% 
% %% Plotting
% 
% 
% figure(1)
% plot(real_times,ent_spec);



