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
        
        function mps_out = apply_to_state_suzuki(obj,mps,chi,time)
            error('Suzuki implementation currently does not work');
            assert(mod(obj.cell_size,2) == 0,'Suzuki integrator only supports even cell size');
            mps_out = mps;
%             times = ones(1,2)./(4 - (4.^(1/3)))*time;
%             times(3) = time - 2*times(1) - 2*times(2); %Suzuki symplectic times
%             times(4) = times(2);
%             times(5) = times(1);
            times = time;
            for j = 1:numel(times)
                mps_out = obj.apply_to_state(mps_out,chi,times(j)/2,[1]);
                mps_out = obj.apply_to_state(mps_out,chi,times(j),[2]);
                mps_out = obj.apply_to_state(mps_out,chi,times(j)/2,[1]);
            end
        end
        
        function varargout = apply_to_state(obj,mps,chi,time,varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            assert(isa(mps,'SpinMPS'));
            assert(obj.dim == mps.dim);
            assert(obj.cell_size == mps.cell_size);
            
            new_states = cell(size(mps.states));
            new_eigs = cell(size(mps.eigs));
            
            apply_op = reshape(expm(-time*obj.log_op),obj.dim,obj.dim,obj.dim,obj.dim);
            curr_state = mps;
            
            if numel(varargin) > 0
                bonds = varargin{1};
            else
                bonds = 1:2;
            end
            
            CHI_TOLERANCE = 1.e-10; %Threshold below which eigs are removed
            max_error = CHI_TOLERANCE;
            
            for k = bonds
                for j = k:2:obj.cell_size
                    ia = j;
                    ib = mod(j,obj.cell_size) + 1;

                    GamA = curr_state.states{ia}; GamB = curr_state.states{ib}; %States on the two sites
                    LamA = curr_state.eigs{ia}; LamB = curr_state.eigs{ib}; %Eigs on the two sites

                    chia = size(GamA,2); chib = size(GamB,3); %Dimensions of the extruding legs
                    
                    assert(size(GamA,3) == size(GamB,2));

                    th = tensorprod(GamA,[1,3,-1],GamB,[2,-1,4]); %Before U operation
                    %th = etprod('abcd',GamA,'ack',GamB,'bkd');
                    th = tensorprod(apply_op,[1,2,-1,-2],th,[-1,-2,3,4]); %Apply the U operator
                    th = tensorprod(diag(LamA),[3,-1],th,[1,2,-1,4]); %Multiply by the old schmidt eigs on the left
                    %At this point, th has 4 indices. First two are vertical
                    %legs, second two are the horizontal extruding legs

                    th = reshape(permute(th,[1,3,2,4]),(obj.dim)*chia,(obj.dim)*chib);
                    [u,S,v] = svd(th); v = v';
                    [evals,ind] = sort(diag(S),'descend');
                    
                    relevant_es = sum(evals > CHI_TOLERANCE);
                    new_chi = min([relevant_es,chi]);
                    if new_chi < relevant_es
                        max_error = max(max_error,evals(ind(new_chi+1)));
                    end
                    new_inds = ind(1:new_chi); %Indices of the eigs we will keep
                    u = u(:,new_inds); v = v(new_inds,:); evals = evals(1:new_chi);
                    evals = evals/sqrt(sum(evals.^2)); %Renormalizing to prevent decay

                    u = reshape(u,obj.dim,chia,new_chi);
                    v = permute(reshape(v,new_chi,obj.dim,chib),[2,1,3]);

                    u = tensorprod(diag(LamA.^(-1)),[-1,2],u,[1,-1,3]); %Multiply by s_a^-1 on left
                    u = tensorprod(diag(evals),[-1,3],u,[1,2,-1]); %Multiply by s_b on the right to restore form

                    new_eigs{ia} = LamA; new_eigs{ib} = evals;
                    new_states{ia} = u; new_states{ib} = v;

                    curr_state = SpinMPS(mps.dim,new_eigs,new_states);
                end
            end
            
            varargout{1} = curr_state;
            if nargout > 1
                varargout{2} = max_error;
            end
            
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
        
        function [bond,hamiltonian] = MPS_Ising_Bond(coupling,field,cell_size)
            dim = 2;
            sz = [[1,0];[0,-1]];
            sx = [[0,1];[1,0]];
            log_op = -coupling*kron(sz,sz) + field*kron(sx,eye(dim));
            bond = MPS_BondOperator(dim,log_op,cell_size);
            hamiltonian = cell(cell_size,1);
            for j = 1:cell_size
                hamiltonian{j} = reshape(log_op,dim,dim,dim,dim);
            end
        end
        
        function [bond,hamiltonian] = MPS_AKLT_Bond(cell_size)
            dim = 3;
            s = cell(1,3);
            s{1} = sqrt(0.5)*[[0,1,0];[1,0,1];[0,1,0]];
            s{2} = sqrt(0.5)*[[0,-1i,0];[1i,0,-1i];[0,1i,0]];
            s{3} = diag([1,0,-1]);
            log_op = zeros(dim^2);
            for j = 1:dim
                log_op = log_op + kron(s{j},s{j});
                for k = 1:dim
                    log_op = log_op + kron(s{j}*s{k},s{j}*s{k})/3;
                end
            end
            bond = MPS_BondOperator(dim,log_op,cell_size);
            hamiltonian = cell(cell_size,1);
            for j = 1:cell_size
                hamiltonian{j} = reshape(log_op,dim,dim,dim,dim);
            end
        end
        
        function [bond,hamiltonian] = MPS_Haldane_Bond(heis,zz,inversion,rotation,xfield,zfield,cell_size)
            dim = 3;
            s = cell(1,3);
            s{1} = sqrt(0.5)*[[0,1,0];[1,0,1];[0,1,0]];
            s{2} = sqrt(0.5)*[[0,-1i,0];[1i,0,-1i];[0,1i,0]];
            s{3} = diag([1,0,-1]);
            log_op = zz*kron(s{3}*s{3},eye(dim)) + xfield*kron(s{1},eye(dim)) + zfield*kron(s{3},eye(dim));
            for j = 1:dim
                log_op = log_op + heis*kron(s{j},s{j});
            end
            log_op = log_op + inversion*kron(s{1}*s{2},s{3}*s{2});
            log_op = log_op + inversion*kron(s{2}*s{1},s{2}*s{3});
            log_op = log_op - inversion*kron(s{3}*s{2},s{1}*s{2});
            log_op = log_op - inversion*kron(s{2}*s{3},s{2}*s{1});
            
            log_op = log_op + rotation*kron(s{1},s{3});
            
            assert(sum(sum(abs(log_op - log_op'))) < 1.e-10,'Hamiltonian must be Hermitian');
            
            e = eye(dim^2);
            inv_op = e([1,4,7,2,5,8,3,6,9],:);
            (inv_op * log_op * inv_op') - log_op
            
            rotz = zeros(dim); rotz(1,3) = 1; rotz(2,2) = 1; rotz(3,1) = 1;
            zop = kron(rotz,rotz);
            (zop * log_op * zop') - log_op
            
            bond = MPS_BondOperator(dim,log_op,cell_size);
            hamiltonian = cell(cell_size,1);
            for j = 1:cell_size
                hamiltonian{j} = reshape(log_op,dim,dim,dim,dim);
            end
        end
    end
end

