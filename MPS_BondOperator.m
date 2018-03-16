classdef MPS_BondOperator
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dim
        log_op
        cell_size
    end
    
    methods
        function obj = MPS_BondOperator(dim,log_op,cell_size)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.dim = dim;
            obj.log_op = log_op;
            obj.cell_size = cell_size;
        end
        
        function mps_out = apply_to_state(obj,mps,chi,time)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            assert(isa(mps,'SpinMPS'));
            assert(obj.dim == mps.dim);
            assert(obj.cell_size == mps.cell_size);
            
            new_states = cell(size(mps.states));
            new_eigs = cell(size(mps.eigs));
            
            apply_op = reshape(expm(time*obj.log_op),obj.dim,obj.dim,obj.dim,obj.dim);
            curr_state = mps;
            
            for k = 1:2
                for j = k:2:obj.cell_size
                    ia = j;
                    ib = mod(j,obj.cell_size) + 1;

                    GamA = curr_state.states{ia}; GamB = curr_state.states{ib}; %States on the two sites
                    LamA = curr_state.eigs{ia}; LamB = curr_state.eigs{ib}; %Eigs on the two sites

                    chia = size(GamA,2); chib = size(GamB,3); %Dimensions of the extruding legs

                    assert(size(GamA,3) == size(GamB,2));

                    th = tensorprod(GamA,[1,3,-1],GamB,[2,-1,4]); %Before U operation
                    %th = etprod('abcd',GamA,'ack',GamB,'bkd');
                    th = tensorprod(th,[-1,-2,3,4],apply_op,[-1,-2,1,2]); %Apply the U operator
                    th = tensorprod(diag(LamA),[3,-1],th,[1,2,-1,4]); %Multiply by the old schmidt eigs on the left
                    %At this point, th has 4 indices. First two are vertical
                    %legs, second two are the horizontal extruding legs

                    th = reshape(permute(th,[1,3,2,4]),(obj.dim)*chia,(obj.dim)*chib);
                    [u,S,v] = svd(th);
                    evals = sort(diag(S));

                    CHI_TOLERANCE = 1.e-10; %Threshold below which eigs are removed
                    new_chi = min(chi,sum(evals > CHI_TOLERANCE));
                    new_inds = 1:new_chi; %Indices of the eigs we will keep
                    u = u(:,new_inds); evals = evals(new_inds); v = v(new_inds,:);
                    evals = evals/sqrt(sum(evals.^2)); %Renormalizing to prevent decay

                    u = reshape(u,obj.dim,chia,new_chi);
                    v = reshape(v,obj.dim,new_chi,chib);

                    u = tensorprod(diag(LamA.^(-1)),[-1,2],u,[1,-1,3]); %Multiply by s_a^-1 on left
                    u = tensorprod(diag(evals),[-1,3],u,[1,2,-1]); %Multiply by s_b on the right to restore form

                    new_eigs{ia} = LamA; new_eigs{ib} = evals;
                    new_states{ia} = u; new_states{ib} = v;

                    curr_state = SpinMPS(mps.dim,new_eigs,new_states);
                end
            end
            
            mps_out = curr_state;
            
        end
    end
    
    methods (Static)
        function out = dot_with_vector(tens,vec,dir)
            vec = reshape(vec,numel(vec));
            ord = 1:numel(size(tens));
            ord(dir) = [];
            ord = [dir,ord];
            tens = permute(tens,ord);
            out = vec .* tens;
            out = permute(out,ord);
        end
        
        function [bond,hamiltonian] = MPS_IsingBond(coupling,field,cell_size)
            dim = 2;
            sz = [[1,0];[0,-1]];
            sx = [[0,1];[1,0]];
            log_op = -coupling*kron(sz,sz) + field*kron(sx,eye(2));
            bond = MPS_BondOperator(dim,log_op,cell_size);
            hamiltonian = cell(cell_size,1);
            for j = 1:cell_size
                hamiltonian{j} = reshape(log_op,dim,dim,dim,dim);
            end
        end
    end
end

