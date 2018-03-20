classdef SpinMPS
    %Matrix product state of a translationally invariant spin chain, with
    %[chi] Schmidt eigenvalues kept.
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        chi
        dim
        eigs
        states %Combination of states AND eigs
        cell_size
    end
    
    methods
        function obj = SpinMPS(dim,eigs,states)
            assert(numel(eigs) == size(eigs,1)); %Column vec
            assert(iscell(eigs) && iscell(states));
            assert(numel(eigs) == numel(states));
            obj.dim = dim;
            obj.eigs = eigs; %Cell of size [cell] with column vectors as components
            obj.states = states; %Cell of size [cell] with 3-tensors as components
            obj.cell_size = numel(eigs);
        end
        
        function vals = bond_expectation(obj,bond_operators)
            assert(iscell(bond_operators));
            assert(numel(bond_operators) == obj.cell_size);
            
            
            vals = cell(size(bond_operators));
            
            for j_bond = 1:numel(bond_operators)
                op = bond_operators{j_bond};
                assert(all(size(op) == ones(1,4)*obj.dim));
                ia = j_bond;
                ib = mod(j_bond,obj.cell_size) + 1;
                
                GamA = obj.states{ia}; GamB = obj.states{ib};
                LamA = obj.eigs{ia};
                
                toprow = tensorprod(GamA,[1,3,-1],GamB,[2,-1,4]); %Multiply the Gammas
                toprow = tensorprod(diag(LamA),[3,-1],toprow,[1,2,-1,4]); %Include the schmidt values on the left
                bottomrow = conj(toprow); %Bottom row - Hermitian conjugate of states
                toprow = tensorprod(toprow,[-1,-2,3,4],op,[-1,-2,1,2]); %Apply the bond operator
                
                vals{j_bond} = trace(tensorprod(toprow,[-1,-2,-3,1],bottomrow,[-1,-2,-3,1]));
                
            end
        end
        
        function spec = entanglement_spectrum(obj,cut_size)
            cells = cut_size/obj.cell_size;
            chil = size(obj.states{1},2); chir = size(obj.states{obj.cell_size},3);
            topcell = obj.states{1};
            celldim = obj.dim^(obj.cell_size);
            for j = 2:obj.cell_size
                topcell = tensorprod(topcell,[1:(j-1),(j+1),-1],obj.states{1},[j,-1,j+2]);
            end
            topcell = reshape(topcell,celldim,chil,chir);
            
            toprow = topcell;
            for j = 2:cells
                toprow = tensorprod(toprow,[1:(j-1),(j+1),-1],topcell,[j,-1,j+2]);
            end
            toprow = tensorprod(diag(obj.eigs{1}),[cells+1,-1],toprow,[1:cells,-1,cells+2]);
            
            densmat = reshape(tensorprod(toprow,[1:cells,-1,-2],conj(toprow),[(cells+1):(2*cells),-1,-2]),celldim.^cells,celldim.^cells);
            spec = eig(densmat);
        end
        
        function test_inversion_symmetry(obj)
            for c = 1:obj.cell_size
                st_c = obj.states{c};
                for j = 1:obj.dim
                    gam = reshape(st_c(j,:,:),size(st_c,2),size(st_c,3));
                    gat = gam.';
                end
            end
        end
    end
    
    methods (Static)
        function mps = initialize_ferromagnet(cell_size,dim)
            eigs = cell(cell_size,1);
            states = cell(cell_size,1);
            
            for j = 1:cell_size
                eigs{j} = [1.0]; %Product state on-site
                sj = zeros(dim,1,1); %The state with non-zero eig
                sj(1,1,1) = 1; sj(2,1,1) = 0; sj(3,1,1) = 0;
                states{j} = sj;
            end
            
            mps = SpinMPS(dim,eigs,states);
        end
    end
end

